
#include "elemat.h"

#include "ssection.h"
#include "material.h"

#include "input/parsers/abaqus/abaqus.h"

using namespace espreso;


Elemat::Elemat()
{

}

Elemat& Elemat::create_dict(const std::vector<SSection> &ssection,const std::vector<AbaqusMaterial> &material)
{
    std::string namee(ssection[0].Elset);
    std::string namee1(ssection[0].Material);

    //for (size_t i=0;i<ssection.size();i++){
    //      std::string name_elset(ssection[i].Elset);
    //      std::string name_mat(ssection[i].Material);

//          elset_mat_dict.emplace(name_elset,name_mat);
    //    elset_mat_dict.insert(std::make_pair(name_elset,name_mat));
//    }
//    int ii=2; int iii=3;
    //for (size_t i=0;i<material.size();i++){
        //std::string name_elset(material[i].Name);
        //mat_int_dict.insert(std::make_pair(ii,iii));
        //ii++;
        //iii++;
        //}
//elset_mat_dict[namee] = namee1;
    return *this;
}
