#ifndef MESHWIDGET_H
#define MESHWIDGET_H

#include "gui/parallel/mpimanager.h"

#include <QtGui>
#include <QOpenGLWidget>

namespace espreso {

    struct MeshRegion
    {
        QVector<float> points;
        QVector3D color;
        bool isActive = true;
    };

    struct MeshMousePosition
    {
        int x;
        int y;
    };

    class MeshWidget : public QOpenGLWidget, protected QOpenGLFunctions
    {
        Q_OBJECT

    public:
        MeshWidget(QWidget* parent = 0);
        MeshWidget(MpiManager* manager, QWidget* parent = 0);
        ~MeshWidget();

        static void initOGL();

        QList<QString> regions();
        void changeRegionState(const QString&);
        QVector3D pickColor(float angle);

    signals:
        void regionClicked(const QString& region);

    protected:
        void initializeGL() override;
        void paintGL() override;
        void resizeGL(int w, int h) override;

        void mouseMoveEvent(QMouseEvent *event) override;
        void wheelEvent(QWheelEvent *event) override;
        void mousePressEvent(QMouseEvent *event) override;

    private:
        MpiManager* m_manager;

        QMap<QString, MeshRegion> m_regions;
        void gatherRegions();

        QOpenGLShaderProgram* m_basicProgram = nullptr;
        GLuint m_basicProgram_position;
        GLuint m_basicProgram_normal;
        QVector3D m_basicProgram_objectColor;
        QVector3D m_basicProgram_lightColor;
        QVector3D m_basicProgram_lightPos;

        float m_lastX;
        float m_lastY;
        float m_viewRotX = 0;
        float m_viewRotY = 0;
        float m_mouse = true;

        static constexpr float FOV = 80.0f;

        float m_fov = FOV;

        QOpenGLShaderProgram* m_clickProgram = nullptr;
        GLuint m_clickProgram_position;
        bool m_clicked = false;
        MeshMousePosition m_mouse_pos;
        void clicked(const MeshMousePosition& event);

        const char* m_basicVS =
                "#version 330 core \n"
                "in highp vec3 aPos;\n"
                "in highp vec3 aNormal;\n"
                "uniform mat4 view;\n"
                "uniform mat4 projection;\n"
                "uniform mat4 model;\n"
                "out vec3 Normal;"
                "out vec3 FragPos;"
                "void main()\n"
                "{\n"
                "gl_Position = projection * view * model * vec4(aPos, 1.0f);\n"
                "FragPos = vec3(gl_Position);\n"
                "Normal = mat3(transpose(inverse(model))) * aNormal;\n"
                "}\n";

        const char* m_basicFS =
                "#version 330 core \n"
                "in vec3 Normal;\n"
                "in vec3 FragPos;\n"
                "uniform vec3 objectColor;\n"
                "uniform vec3 lightColor;\n"
                "uniform vec3 lightPos;\n"
                "out vec4 FragColor;\n"
                "void main()\n"
                "{\n"
                "float ambientStrength = 0.2;\n"
                "vec3 ambient = ambientStrength * lightColor;\n"
                "vec3 norm = normalize(Normal);\n"
                "vec3 lightDir = normalize(lightPos - FragPos);\n"
                "float diff = max(dot(norm, lightDir), 0.0);\n"
                "vec3 diffuse = diff * lightColor;\n"
                "vec3 result = (ambient + diffuse) * objectColor;\n"
                "FragColor = vec4(result, 1.0);\n"
                "}\n";

        const char* m_clickVS =
                "#version 330 core\n"
                "uniform mat4 projection;\n"
                "uniform mat4 view;\n"
                "in vec4 position;\n"
                "void main()\n"
                "{\n"
                "gl_Position = projection * view * position;"
                "}\n";
        const char* m_clickFS =
                "#version 330\n"
                "uniform vec3 code;\n"
                "out vec4 outputF;\n"
                "void main()\n"
                "{\n"
                "outputF = vec4(code.x / 255.0f, code.y / 255.0f, code.z / 255.0f, 0);"
                "}\n";
    };

}

#endif // MESHWIDGET_H
