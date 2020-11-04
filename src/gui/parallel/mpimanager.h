
#ifndef MPIMANAGER_H
#define MPIMANAGER_H

#include <QtGui>

namespace espreso
{

struct ECF;
class Mesh;


namespace gui
{
    class MpiOperation
    {
    public:
        static const int EXIT = 0;
        static const int OPEN_ECF = 1;
        static const int GATHER_MESH = 2;
    };
}

class MpiManager
{
public:
    MpiManager(int argc, char* argv[]);
    ~MpiManager();
    ECF* ecf();

    Mesh* mesh();

    void loop();

    void masterOpenECF(const QString&);
    QMap<QString, QVector<float> >* masterGatherMesh();
    void masterExit();


private:
    float m_threshold = 5.0f;

    void performOperation(int code);

    void slaveOpenECF();
    void slaveGatherMesh();

    void _openECF(const std::string& filename);
    QMap<QString, QVector<float> >* _gatherMesh();
    QVector3D pickColor(float angle);
};

}

#endif // MPIMANAGER_H
