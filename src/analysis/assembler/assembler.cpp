
#include "assembler.hpp"

#include "basis/evaluator/evaluator.h"
#include "basis/utilities/parser.h"
#include "basis/utilities/utils.h"
#include "config/holders/expression.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/eslog.hpp"
#include "esinfo/meshinfo.h"
#include "mesh/store/domainstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include "wrappers/mpi/communication.h"
#include "wrappers/exprtk/exprtk.h"

#include <algorithm>
#include <numeric>

#include <iostream>

using namespace espreso;

Assembler::Assembler(PhysicsConfiguration &settings)
: settings(settings), withBEM(false), threaded(false)
{
    BEM.resize(info::mesh->domains->size);
    BETI.resize(info::mesh->domains->size);
    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
            for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
                info::mesh->elements->eintervals[i].thread = t;
            }
        }
    }
    for (int t = 0; t < info::env::threads; ++t) {
        for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
            for (size_t r = 0; r < info::mesh->boundary.size(); ++r) {
                if (info::mesh->boundary[r]->dimension) {
                    for (esint i = info::mesh->boundary[r]->eintervalsDistribution[d]; i < info::mesh->boundary[r]->eintervalsDistribution[d + 1]; ++i) {
                        info::mesh->boundary[r]->eintervals[i].thread = t;
                    }
                }
            }
        }
    }
}

Assembler::~Assembler()
{

}

void Assembler::assemble(const SubKernel::Action action)
{
    if (threaded) {
        #pragma omp parallel for
        for (int t = 0; t < info::env::threads; ++t) {
            for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
                if (BEM[d] && action == SubKernel::Action::ASSEMBLE) {
                    bem(action, d, BETI[d]);
                } else {
                    for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
                        elements(action, i);
                    }
                }
            }
        }
        #pragma omp parallel for
        for (int t = 0; t < info::env::threads; ++t) {
            for (size_t r = 1; r < info::mesh->boundary.size(); ++r) {
                if (info::mesh->boundary[r]->dimension) {
                    for (size_t d = info::mesh->domains->distribution[t]; d < info::mesh->domains->distribution[t + 1]; d++) {
                        for (esint i = info::mesh->boundary[r]->eintervalsDistribution[d]; i < info::mesh->boundary[r]->eintervalsDistribution[d + 1]; ++i) {
                            boundary(action, r, i);
                        }
                    }
                }
            }
        }
        for (size_t r = 1; r < info::mesh->boundary.size(); ++r) {
            for (int t = 0; t < info::env::threads; ++t) {
                nodes(action, r, t); // never parallel
            }
        }

    } else {
        for (esint d = 0; d < info::mesh->domains->size; d++) {
            if (BEM[d] && action == SubKernel::Action::ASSEMBLE) {
                bem(action, d, BETI[d]);
            } else {
                for (esint i = info::mesh->elements->eintervalsDistribution[d]; i < info::mesh->elements->eintervalsDistribution[d + 1]; ++i) {
                    elements(action, i);
                }
            }
        }
        for (size_t r = 1; r < info::mesh->boundary.size(); ++r) {
            if (info::mesh->boundary[r]->dimension) {
                for (esint d = 0; d < info::mesh->domains->size; d++) {
                    for (esint i = info::mesh->boundary[r]->eintervalsDistribution[d]; i < info::mesh->boundary[r]->eintervalsDistribution[d + 1]; ++i) {
                        boundary(action, r, i);
                    }
                }
            }
            for (int t = 0; t < info::env::threads; ++t) {
                nodes(action, r, t);
            }
        }
    }
}

void Assembler::printElementVolume(std::vector<double> &volume)
{
    std::vector<double> sum(info::mesh->elementsRegions.size());
    for (size_t r = 1; r < info::mesh->elementsRegions.size(); ++r) {
        for (size_t i = 0; i < info::mesh->elements->eintervals.size(); ++i) {
            if (info::mesh->elements->eintervals[i].region == (esint)r || (info::mesh->elements->eintervals[i].region == 0 && r == info::mesh->elementsRegions.size() - 1)) {
                sum[0] += volume[i];
                sum[r] += volume[i];
            }
        }
    }
    Communication::allReduce(sum, Communication::OP::SUM);

    eslog::info("  ELEMENT REGION VOLUME                                                                        \n");
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
    for (size_t r = 0; r < info::mesh->elementsRegions.size(); ++r) {
        eslog::info("     %30s :                                            %e   \n", info::mesh->elementsRegions[r]->name.c_str(), sum[r]);
    }
}

void Assembler::printBoundarySurface(std::vector<double> &surface)
{
    Communication::allReduce(surface, Communication::OP::SUM);

    eslog::info("\n  BOUDNARY REGION SURFACE                                                                      \n");
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - \n");
    for (size_t r = 1; r < info::mesh->boundaryRegions.size(); ++r) {
        if (info::mesh->boundaryRegions[r]->dimension) {
            eslog::info("     %30s :                                            %e   \n", info::mesh->boundaryRegions[r]->name.c_str(), surface[r]);
            info::mesh->boundaryRegions[r]->area = surface[r];
        }
    }
}

void Assembler::printMaterials(const std::map<std::string, std::string> &settings)
{
    for (auto reg = info::mesh->elementsRegions.begin() + 1; reg != info::mesh->elementsRegions.end(); ++reg) {
        auto ms = settings.find((*reg)->name);
        if (ms != settings.end()) {
            eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
        } else {
            ms = settings.find("ALL_ELEMENTS");
            if (ms != settings.end()) {
                eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), ms->second.c_str());
            } else {
                eslog::info("  %55s: %34s\n", (*reg)->name.c_str(), info::mesh->materials.front()->name.c_str());
            }
        }
    }
}

bool Assembler::checkExpression(const std::string &name, ECFExpression &expression)
{
    if (expression.evaluator == nullptr) {
        if (Exprtk::check(expression.value)) {
            eslog::warning("   %25s:  %62s \n", name.c_str(), "INVALID EXPRESSION");
            return false;
        }
        expression.evaluator = Evaluator::create(expression.value);
    }
    if (expression.evaluator->parameters.size()) {
        eslog::info("   %25s:  %62s \n", name.c_str(), expression.value.c_str());
    } else {
        eslog::info("   %25s:  %62g \n", name.c_str(), expression.evaluator->evaluate());
    }
    return true;
}

bool Assembler::checkExpression(const std::string &name, ECFExpressionVector &expression)
{
    bool correct = true;
    correct &= checkExpression(name + ".X", expression.x);
    correct &= checkExpression(name + ".Y", expression.y);
    if (info::mesh->dimension == 3) {
        correct &= checkExpression(name + ".Z", expression.z);
    }
    return correct;
}

bool Assembler::checkElementParameter(const std::string &name, std::map<std::string, ECFExpression> &settings)
{
    for (auto s = settings.begin(); s != settings.end(); ++s) {
        info::mesh->bregion(s->first);
    }
    if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
        return checkExpression(name, settings.begin()->second);
    } else {
        eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
        for (auto region = info::mesh->elementsRegions.crbegin(); region != info::mesh->elementsRegions.crend(); ++region) {
            auto it = settings.find((*region)->name);
            if (it != settings.end()) {
                if (!checkExpression(it->first, it->second)) { return false; }
            }
        }
        eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    }
    return true;
}

bool Assembler::checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings)
{
    for (auto s = settings.begin(); s != settings.end(); ++s) {
        info::mesh->bregion(s->first);
    }
    if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
        switch (info::mesh->dimension) {
        case 2: if (!checkExpression(name + ".X", settings.begin()->second.x) || !checkExpression(name + ".Y", settings.begin()->second.y)) { return false; } break;
        case 3: if (!checkExpression(name + ".X", settings.begin()->second.x) || !checkExpression(name + ".Y", settings.begin()->second.y) || !checkExpression(name + ".Z", settings.begin()->second.z)) { return false; } break;
        }
    } else {
        eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
        for (auto region = info::mesh->elementsRegions.crbegin(); region != info::mesh->elementsRegions.crend(); ++region) {
            auto it = settings.find((*region)->name);
            if (it != settings.end()) {
                switch (info::mesh->dimension) {
                case 2: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y)) { return false; } break;
                case 3: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y) || !checkExpression(it->first + ".Z", it->second.z)) { return false; } break;
                }
            }
        }
    }
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    return true;
}

bool Assembler::checkElementParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings, int dim)
{
    for (auto s = settings.begin(); s != settings.end(); ++s) {
        info::mesh->bregion(s->first);
    }
    if (settings.size() == 1 && StringCompare::caseInsensitiveEq(settings.begin()->first, "ALL_ELEMENTS")) {
        return checkExpression(name, settings.begin()->second.data[dim]);
    } else {
        eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
        for (auto region = info::mesh->elementsRegions.crbegin(); region != info::mesh->elementsRegions.crend(); ++region) {
            auto it = settings.find((*region)->name);
            if (it != settings.end()) {
                return checkExpression(name, it->second.data[dim]);
            }
        }
        eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    }
    return true;
}

bool Assembler::checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpression> &settings)
{
    for (auto s = settings.begin(); s != settings.end(); ++s) {
        info::mesh->bregion(s->first);
    }
    eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
    for (auto region = info::mesh->boundaryRegions.crbegin(); region != info::mesh->boundaryRegions.crend(); ++region) {
        auto it = settings.find((*region)->name);
        if (it != settings.end()) {
            if (!checkExpression(it->first, it->second)) { return false; }
        }
    }
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    return true;
}

bool Assembler::checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionVector> &settings)
{
    for (auto s = settings.begin(); s != settings.end(); ++s) {
        info::mesh->bregion(s->first);
    }
    eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
    for (auto region = info::mesh->boundaryRegions.crbegin(); region != info::mesh->boundaryRegions.crend(); ++region) {
        auto it = settings.find((*region)->name);
        if (it != settings.end()) {
            switch (info::mesh->dimension) {
            case 2: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y)) { return false; } break;
            case 3: if (!checkExpression(it->first + ".X", it->second.x) || !checkExpression(it->first + ".Y", it->second.y) || !checkExpression(it->first + ".Z", it->second.z)) { return false; } break;
            }
        }
    }
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    return true;
}

bool Assembler::checkBoundaryParameter(const std::string &name, std::map<std::string, ECFExpressionOptionalVector> &settings)
{
    for (auto s = settings.begin(); s != settings.end(); ++s) {
        info::mesh->bregion(s->first);
    }
    eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
    for (auto region = info::mesh->boundaryRegions.crbegin(); region != info::mesh->boundaryRegions.crend(); ++region) {
        auto it = settings.find((*region)->name);
        if (it != settings.end()) {
            switch (info::mesh->dimension) {
            case 2: if ((it->second.x.isset && !checkExpression(it->first + ".X", it->second.x)) || (it->second.y.isset && !checkExpression(it->first + ".Y", it->second.y))) { return false; } break;
            case 3: if ((it->second.x.isset && !checkExpression(it->first + ".X", it->second.x)) || (it->second.y.isset && !checkExpression(it->first + ".Y", it->second.y)) || (it->second.z.isset && !checkExpression(it->first + ".Z", it->second.z))) { return false; } break;
            }
        }
    }
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    return true;
}

bool Assembler::checkBoundaryParameter(const std::string &name, std::map<std::string, ECFHarmonicExpressionVector> &settings)
{
    for (auto s = settings.begin(); s != settings.end(); ++s) {
        info::mesh->bregion(s->first);
    }
    eslog::info("  %s%*s \n", name.c_str(), 91 - name.size(), "");
    for (auto region = info::mesh->boundaryRegions.crbegin(); region != info::mesh->boundaryRegions.crend(); ++region) {
        auto it = settings.find((*region)->name);
        if (it != settings.end()) {
            switch (info::mesh->dimension) {
            case 2:
                if (
                        !checkExpression(it->first + ".X.MAGNITUDE", it->second.magnitude.x) || !checkExpression(it->first + ".X.PHASE    ", it->second.phase.x) ||
                        !checkExpression(it->first + ".Y.MAGNITUDE", it->second.magnitude.y) || !checkExpression(it->first + ".Y.PHASE    ", it->second.phase.y))
                { return false; }
                break;
            case 3:
                if (
                        !checkExpression(it->first + ".X.MAGNITUDE", it->second.magnitude.x) || !checkExpression(it->first + ".X.PHASE    ", it->second.phase.x) ||
                        !checkExpression(it->first + ".Y.MAGNITUDE", it->second.magnitude.y) || !checkExpression(it->first + ".Y.PHASE    ", it->second.phase.y) ||
                        !checkExpression(it->first + ".Z.MAGNITUDE", it->second.magnitude.z) || !checkExpression(it->first + ".Z.PHASE    ", it->second.phase.z))
                { return false; }
                break;
            }
        }
    }
    eslog::info("  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  \n");
    return true;
}

ECFExpression* Assembler::getExpression(size_t interval, std::map<std::string, ECFExpression> &settings)
{
    int region = info::mesh->elements->eintervals[interval].region;
    auto it = settings.find(info::mesh->elementsRegions[region]->name);
    if (it == settings.end()) {
        it = settings.find(info::mesh->elementsRegions[0]->name);
    }
    return it != settings.end() ? &it->second : nullptr;
}

ECFExpressionVector* Assembler::getExpression(size_t interval, std::map<std::string, ECFExpressionVector> &settings)
{
    int region = info::mesh->elements->eintervals[interval].region;
    auto it = settings.find(info::mesh->elementsRegions[region]->name);
    if (it == settings.end()) {
        it = settings.find(info::mesh->elementsRegions[0]->name);
    }
    return it != settings.end() ? &it->second : nullptr;
}

ECFHarmonicExpressionVector* Assembler::getExpression(size_t interval, std::map<std::string, ECFHarmonicExpressionVector> &settings)
{
    int region = info::mesh->elements->eintervals[interval].region;
    auto it = settings.find(info::mesh->elementsRegions[region]->name);
    if (it == settings.end()) {
        it = settings.find(info::mesh->elementsRegions[0]->name);
    }
    return it != settings.end() ? &it->second : nullptr;
}

ECFExpression* Assembler::getExpression(const std::string &name, std::map<std::string, ECFExpression> &settings)
{
    auto it = settings.find(name);
    return it != settings.end() ? &it->second : nullptr;
}

ECFExpressionVector* Assembler::getExpression(const std::string &name, std::map<std::string, ECFExpressionVector> &settings)
{
    auto it = settings.find(name);
    return it != settings.end() ? &it->second : nullptr;
}

ECFHarmonicExpressionVector* Assembler::getExpression(const std::string &name, std::map<std::string, ECFHarmonicExpressionVector> &settings)
{
    auto it = settings.find(name);
    return it != settings.end() ? &it->second : nullptr;
}

Evaluator* Assembler::getEvaluator(size_t interval, std::map<std::string, ECFExpression> &settings)
{
    int region = info::mesh->elements->eintervals[interval].region;
    auto it = settings.find(info::mesh->elementsRegions[region]->name);
    if (it == settings.end()) {
        it = settings.find(info::mesh->elementsRegions[0]->name);
    }
    if (it != settings.end()) {
        return it->second.evaluator;
    }
    return nullptr;
}

Evaluator* Assembler::getEvaluator(size_t interval, std::map<std::string, ECFExpressionVector> &settings, int dim)
{
    int region = info::mesh->elements->eintervals[interval].region;
    auto it = settings.find(info::mesh->elementsRegions[region]->name);
    if (it == settings.end()) {
        it = settings.find(info::mesh->elementsRegions[0]->name);
    }
    if (it != settings.end()) {
        return it->second.data[dim].evaluator;
    }
    return nullptr;
}

bool Assembler::isBEM(size_t interval)
{
    int region = info::mesh->elements->eintervals[interval].region;
    auto it = settings.discretization.find(info::mesh->elementsRegions[region]->name);
    if (it == settings.discretization.end()) {
        it = settings.discretization.find(info::mesh->elementsRegions[0]->name);
    }
    if (it != settings.discretization.end() && it->second == PhysicsConfiguration::DISCRETIZATION::BEM) {
        return true;
    }
    return false;
}

