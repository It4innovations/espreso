
#include "mpimanager.h"

#include "basis/containers/serializededata.h"
#include "wrappers/mpi/communication.h"
#include "config/ecf/ecf.h"
#include "esinfo/mpiinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "mesh/store/elementstore.h"
#include "mesh/store/nodestore.h"
#include "mesh/store/surfacestore.h"
#include "mesh/store/elementsregionstore.h"
#include "mesh/store/boundaryregionstore.h"

#include <QDebug>
#include <QtGui>
#include <algorithm>

using namespace espreso;

MpiManager::MpiManager(int argc, char* argv[])
{

}

MpiManager::~MpiManager()
{

}

ECF* MpiManager::ecf()
{
    return info::ecf;
}

Mesh* MpiManager::mesh()
{
    return info::mesh;
}

void MpiManager::loop()
{
    if (info::mpi::rank == 0) return;

    while (true)
    {
        int opcode = 0;
        MPI_Bcast(&opcode, 1, MPI_INT, 0, info::mpi::comm);
        if (opcode == 0) break;

        this->performOperation(opcode);
    }
}

void MpiManager::performOperation(int code)
{
    switch (code)
    {
        case gui::MpiOperation::OPEN_ECF:
            this->slaveOpenECF();
            break;
        case gui::MpiOperation::GATHER_MESH:
            this->slaveGatherMesh();
            break;
        default:
            ;
    }
}

void MpiManager::masterOpenECF(const QString& filename)
{
    int op = gui::MpiOperation::OPEN_ECF;
    MPI_Bcast(&op, 1, MPI_INT, 0, info::mpi::comm);

    std::string fname = filename.toStdString();
    std::vector<char> fn(fname.begin(), fname.end());
    fn.push_back(0);
    int len = fname.size();
    MPI_Bcast(&len, 1, MPI_INT, 0, info::mpi::comm);
    MPI_Bcast(&fn[0], len, MPI_CHAR, 0, info::mpi::comm);

    this->_openECF(filename.toStdString());
}

void MpiManager::slaveOpenECF()
{
    int len = 0;
    MPI_Bcast(&len, 1, MPI_INT, 0, info::mpi::comm);
    char fn[len + 1];
    MPI_Bcast(fn, len, MPI_CHAR, 0, info::mpi::comm);
    fn[len] = 0;
    this->_openECF(std::string(fn));
}

void MpiManager::_openECF(const std::string& filename)
{
    ECF::init(filename);
}

QMap<QString, QVector<float> >* MpiManager::masterGatherMesh()
{
    int op = gui::MpiOperation::GATHER_MESH;
    MPI_Bcast(&op, 1, MPI_INT, 0, info::mpi::comm);

    return this->_gatherMesh();
}

void MpiManager::slaveGatherMesh()
{
    this->_gatherMesh();
}

QMap<QString, QVector<float> >* MpiManager::_gatherMesh()
{
	Mesh::init();
	Mesh::load();
	info::mesh->preprocessForGUI();

	QMap<QString, QVector<float> > regions;

    float axis_min = -1.0f;
    float axis_max = 1.0f;
    float axis_len = axis_max - axis_min;

    Point rmin(1e10, 1e10, 1e10), rmax(-1e10, -1e10, -1e10), min, max;

    for (esint i = 0; i < info::mesh->nodes->size; i++) {
        rmin.x = std::min(rmin.x, info::mesh->nodes->coordinates->datatarray()[i].x);
        rmin.y = std::min(rmin.y, info::mesh->nodes->coordinates->datatarray()[i].y);
        rmin.z = std::min(rmin.z, info::mesh->nodes->coordinates->datatarray()[i].z);
        rmax.x = std::max(rmax.x, info::mesh->nodes->coordinates->datatarray()[i].x);
        rmax.y = std::max(rmax.y, info::mesh->nodes->coordinates->datatarray()[i].y);
        rmax.z = std::max(rmax.z, info::mesh->nodes->coordinates->datatarray()[i].z);
    }

    Communication::allReduce(&rmin.x, &min.x, 3, MPI_DOUBLE, MPI_MIN);
    Communication::allReduce(&rmax.x, &max.x, 3, MPI_DOUBLE, MPI_MAX);

    float x_len = max.x - min.x;
    float y_len = max.y - min.y;
    float z_len = max.z - min.z;
    QVector<float> axises;
    axises << x_len << y_len << z_len;
    float max_axis_len = *std::max_element(axises.begin(), axises.end());

    auto transform = [] (float coordinate, float old_axis_min, float old_axis_len, float max_axis,
            float new_axis_len, float new_axis_min)
    {
        return ( (coordinate - old_axis_min + ( (max_axis - old_axis_len) / 2.0f ) ) / max_axis ) * new_axis_len + new_axis_min;
    };

    auto createRegion = [&] (QString name, serializededata<esint, esint>* triangles)
    {

        QVector<float> &mesh = regions[name];

        for (
                auto t = triangles->cbegin();
                t != triangles->cend();
                ++t) {

            for (int i = 0; i < 3; ++i) {
                Point p = info::mesh->nodes->coordinates->datatarray()[t->at(i)];
                float _x = transform(p.x, min.x, x_len, max_axis_len, axis_len, axis_min);
                float _y = transform(p.y, min.y, y_len, max_axis_len, axis_len, axis_min);
                float _z = transform(p.z, min.z, z_len, max_axis_len, axis_len, axis_min);
                mesh.push_back( _x );
                mesh.push_back( _y );
                mesh.push_back( _z );
                mesh.push_back(0);
                mesh.push_back(0);
                mesh.push_back(0);
            }
        }

        struct vertex
        {
            QVector<float*> locations;
            QVector<QVector3D> triangle_normals;
        };

        QHash<QString, vertex> vertices;

        // COMPUTE TRIANGLE NORMALS
        for (int t = 0; t < mesh.size(); t += 18)
        {
            QVector3D v1(mesh[t], mesh[t + 1], mesh[t + 2]);
            QVector3D v2(mesh[t + 6], mesh[t + 7], mesh[t + 8]);
            QVector3D v3(mesh[t + 12], mesh[t + 13], mesh[t + 14]);

            QVector3D edge1 = v2 - v1;
            QVector3D edge2 = v3 - v1;

            QVector3D normal = QVector3D::crossProduct(edge1, edge2).normalized();

            for (int v = 0; v < 3; v++)
            {
                mesh[t + 3 + 6*v] = normal.x();
                mesh[t + 3 + 6*v + 1] = normal.y();
                mesh[t + 3 + 6*v + 2] = normal.z();
            }

            QVector<QVector3D> vs;
            vs << v1 << v2 << v3;

            for (int i = 0; i < vs.size(); i++)
            {
                QString v_key;
                QDebug(&v_key) << v1;
                if (vertices.find(v_key) == vertices.end())
                {
                    struct vertex v;
                    vertices[v_key] = v;
                }
                vertices[v_key].locations.append(&mesh.data()[t + 3 + 6*i]);
                vertices[v_key].triangle_normals.append(normal);
            }
        }

        // COMPUTE VERTEX NORMAL FROM SURROUNDING TRIANGLES AND USE IT WHEN angle IS BELOW m_threshold
        QHashIterator<QString, vertex> i(vertices);
        while (i.hasNext()) {
            i.next();
            QVector3D normal;
            foreach (QVector3D n, i.value().triangle_normals) {
                normal += n;
            }
            normal.normalize();


            for (int v = 0; v < i.value().locations.size(); v++)
            {
                float cos_angle = QVector3D::dotProduct(normal, i.value().triangle_normals[v]) /
                        ( normal.length() * i.value().triangle_normals[v].length() );
                float angle = qRadiansToDegrees(qAcos(cos_angle));
                if (angle < m_threshold)
                {
                    i.value().locations[v][0] = normal.x();
                    i.value().locations[v][1] = normal.y();
                    i.value().locations[v][2] = normal.z();
                }
            }
        }

    };

    for (auto region = info::mesh->elementsRegions.begin() + 1;
         region != info::mesh->elementsRegions.end();
         region++)
    {
        createRegion(QString::fromStdString((*region)->name),
                     (*region)->surface->triangles);

    }

    for (auto region = info::mesh->boundaryRegions.begin();
         region != info::mesh->boundaryRegions.end();
         region++)
    {
        if ((*region)->dimension < 2) continue;
        createRegion(QString::fromStdString((*region)->name),
                     (*region)->triangles);

    }

//    for (size_t e = 0; e < info::mesh->elements().size(); e++) {
//
//        for (size_t f = 0; f < info::mesh->elements()[e]->faces(); f++) {
//            QVector<QString> regionNames;
//            regionNames << QLatin1String("#global");
//            if (info::mesh->elements()[e]->face(f)->regions().size())
//            {
//                regionNames.clear();
//
//                for (size_t r = 0; r < info::mesh->elements()[e]->face(f)->regions().size(); r++)
//                {
//                    QString regionName = QString::fromStdString(info::mesh->elements()[e]->face(f)->regions()[0]->name);
//                    regionNames << regionName;
//
//                    if (!regions.contains(regionName))
//                    {
//                        regions.insert(regionName, QVector<float>());
//                    }
//                }
//            }
//
//            std::vector<std::vector<esint> > triangles = dynamic_cast<PlaneElement*>(info::mesh->elements()[e]->face(f))->triangularize();
//
//            for (size_t t = 0; t < triangles.size(); t++) {
//
//                for (size_t n = 0; n < triangles[t].size(); n++) {
//                    foreach (QString rn, regionNames)
//                    {
//                        regions[rn].push_back(info::mesh->coordinates()[triangles[t][n]].x);
//                        regions[rn].push_back(info::mesh->coordinates()[triangles[t][n]].y);
//                        regions[rn].push_back(info::mesh->coordinates()[triangles[t][n]].z);
//                        regions[rn].push_back(0.0f);
//                        regions[rn].push_back(1.0f);
//                        regions[rn].push_back(0.0f);
//                    }
//                }
//
//            }
//        }
//    }

    QMapIterator<QString, QVector<float> > it(regions);

    QMap<QString, QVector<float> >* r_regions = (info::mpi::rank == 0) ? new QMap<QString, QVector<float> >() : nullptr;

    while (it.hasNext())
    {
        it.next();

        int num = it.value().size();
        int nums[info::mpi::size];
        MPI_Gather(&num, 1, MPI_INT, nums, 1, MPI_INT, 0, info::mpi::comm);

        int numsum = 0;
        int displs[info::mpi::size];
        QVector<float> coordinates;
        if (info::mpi::rank == 0)
        {
            for (int i = 0; i < info::mpi::size; i++)
            {
                displs[i] = numsum;
                numsum += nums[i];
            }
            coordinates.resize(numsum);
        }

        MPI_Gatherv(const_cast<float*>(it.value().data()), num, MPI_FLOAT, coordinates.data(), nums, displs, MPI_FLOAT, 0, info::mpi::comm);

        if (info::mpi::rank == 0)
        {
            r_regions->insert(it.key(), coordinates);
        }
    }

    return r_regions;
}

void MpiManager::masterExit()
{
    int code = gui::MpiOperation::EXIT;
    MPI_Bcast(&code, 1, MPI_INT, 0, info::mpi::comm);
}
