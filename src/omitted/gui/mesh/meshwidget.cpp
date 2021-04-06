#include "meshwidget.h"

#include "config/ecf/environment.h"
#include "mpi.h"

#include <QDebug>
#include <QVBoxLayout>

using namespace espreso;

void MeshWidget::initOGL()
{
    QSurfaceFormat format;
    format.setVersion(3,3);
    format.setProfile(QSurfaceFormat::CoreProfile);
    QSurfaceFormat::setDefaultFormat(format);
}

MeshWidget::MeshWidget(QWidget *parent) : QOpenGLWidget(parent)
{
}

MeshWidget::MeshWidget(MpiManager* manager, QWidget* parent) :
    MeshWidget(parent)
{
    this->m_manager = manager;
    this->gatherRegions();
    this->setLayout(new QVBoxLayout);

//    this->computeMesh();
}

MeshWidget::~MeshWidget()
{
    if (this->m_basicProgram != nullptr)
        delete this->m_basicProgram;
    if (this->m_clickProgram != nullptr)
        delete this->m_clickProgram;
}

void MeshWidget::gatherRegions()
{
    auto regions = this->m_manager->masterGatherMesh();
    QMapIterator<QString, QVector<float> > it(*regions);

    float regionsNum = regions->size();
    float offset = 360.0f / regionsNum;
    float step = 0.0f;

    while (it.hasNext())
    {
        it.next();
        MeshRegion reg;
        reg.color = this->pickColor(step);
        step += offset;
        reg.points = it.value();
        this->m_regions.insert(it.key(), reg);
    }

    delete regions;
}

void MeshWidget::initializeGL()
{
    this->initializeOpenGLFunctions();

    this->m_basicProgram = new QOpenGLShaderProgram(this);
    m_basicProgram->addShaderFromSourceCode(QOpenGLShader::Vertex, m_basicVS);
    m_basicProgram->addShaderFromSourceCode(QOpenGLShader::Fragment, m_basicFS);
    m_basicProgram->link();

    this->m_basicProgram_position = m_basicProgram->attributeLocation("aPos");
    this->m_basicProgram_normal = m_basicProgram->attributeLocation("aNormal");

    this->m_basicProgram_objectColor = QVector3D(1.0f, 0.5f, 0.31f);
    this->m_basicProgram_lightColor = QVector3D(1.0f, 1.0f, 1.0f);
    this->m_basicProgram_lightPos = QVector3D(0.0f, 0.0f, 5.0f);

    this->m_lastX = width() / 2;
    this->m_lastY = height() / 2;

    this->m_clickProgram = new QOpenGLShaderProgram(this);
    m_clickProgram->addShaderFromSourceCode(QOpenGLShader::Vertex, m_clickVS);
    m_clickProgram->addShaderFromSourceCode(QOpenGLShader::Fragment, m_clickFS);
    m_clickProgram->link();

    this->m_clickProgram_position = m_clickProgram->attributeLocation("position");

    glEnable(GL_DEPTH_TEST);
}

void MeshWidget::paintGL()
{
    if (this->m_clicked)
    {
        this->m_clicked = false;
        this->clicked(this->m_mouse_pos);
    }

    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_basicProgram->bind();
    m_basicProgram->setUniformValue("objectColor", m_basicProgram_objectColor);
    m_basicProgram->setUniformValue("lightColor", m_basicProgram_lightColor);
    m_basicProgram->setUniformValue("lightPos", m_basicProgram_lightPos);

    QMatrix4x4 view;
    view.translate(0.0f, 0.0f, -3.0f);
    m_basicProgram->setUniformValue("view", view);

    QMatrix4x4 projection;
    projection.perspective(m_fov, (float)width() / (float)height(), 0.1f, 100.0f);
    m_basicProgram->setUniformValue("projection", projection);

    foreach (MeshRegion r, m_regions) {

        if (!r.isActive) continue;

        m_basicProgram->setUniformValue("objectColor", r.color);

        float* vertices = &r.points.data()[0];
        int len = r.points.size() / 6;

        QMatrix4x4 model;
        model.rotate(m_viewRotX, 1.0f, 0.0f, 0.0f);
        model.rotate(m_viewRotY, 0.0f, 1.0f, 0.0f);
        this->m_basicProgram->setUniformValue("model", model);

        QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
        vbo.create();
        vbo.bind();
        vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
        vbo.allocate(vertices, sizeof(float) * r.points.size());

        QOpenGLVertexArrayObject cubeVAO(this);
        cubeVAO.create();
        cubeVAO.bind();

        glVertexAttribPointer(m_basicProgram_position, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glVertexAttribPointer(m_basicProgram_normal, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(1);

//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glDrawArrays(GL_TRIANGLES, 0, len);

        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);

        cubeVAO.release();

        vbo.release();
    }

    m_basicProgram->release();
}

void MeshWidget::resizeGL(int w, int h)
{

}

void MeshWidget::mouseMoveEvent(QMouseEvent *event)
{
    float xpos = event->pos().x();
    float ypos = event->pos().y();

    if (m_mouse)
    {
        this->m_lastX = xpos;
        this->m_lastY = ypos;
        m_mouse = false;
    }

    float xoffset = xpos - m_lastX;
    float yoffset = m_lastY - ypos;
    this->m_lastX = xpos;
    this->m_lastY = ypos;

    float sensitivity = 0.35f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    this->m_viewRotX -= yoffset;
    this->m_viewRotY += xoffset;

    if (m_viewRotX > 179.0f)
        m_viewRotX = 179.0f;
    if (m_viewRotX < -179.0f)
        m_viewRotX = -179.0f;

    if (m_viewRotY > 179.0f)
        m_viewRotY = 179.0f;
    if (m_viewRotY < -179.0f)
        m_viewRotY = -179.0f;

    this->update();
}

void MeshWidget::wheelEvent(QWheelEvent *event)
{
    float sign = (event->angleDelta().y() > 0) ? 1.0f : -1.0f;

    if (m_fov >= 1.0f && m_fov <= FOV)
        m_fov -= sign * 3;

    if (m_fov <= 1.0f)
        m_fov = 1.0f;

    if (m_fov >= FOV)
        m_fov = FOV;

    this->update();
}

void MeshWidget::clicked(const MeshMousePosition& position)
{
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    m_clickProgram->bind();

    QMatrix4x4 view;
    view.translate(0.0f, 0.0f, -3.0f);
    view.rotate(m_viewRotX, 1.0f, 0.0f, 0.0f);
    view.rotate(m_viewRotY, 0.0f, 1.0f, 0.0f);

    m_clickProgram->setUniformValue("view", view);

    QMatrix4x4 projection;
    projection.perspective(m_fov, (float)width() / (float)height(), 0.1f, 100.0f);

    m_clickProgram->setUniformValue("projection", projection);

    int color = 1;
    QHash<int, QString> colorCodes;

    QMapIterator<QString, MeshRegion> it(m_regions);

    while (it.hasNext()) {
        it.next();
        if (!it.value().isActive) continue;

        MeshRegion r = it.value();

        float red = (color & 0x00FF0000) >> 16;
        float green = (color & 0x0000FF00) >> 8;
        float blue = (color & 0x000000FF);
        m_clickProgram->setUniformValue("code", QVector3D(red, green, blue));
        colorCodes.insert(color, it.key());
        color += 1;

        float* vertices = &r.points.data()[0];
        int len = r.points.size() / 6;

        QOpenGLBuffer vbo(QOpenGLBuffer::VertexBuffer);
        vbo.create();
        vbo.bind();
        vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
        vbo.allocate(vertices, sizeof(float) * r.points.size());

        QOpenGLVertexArrayObject vao(this);
        vao.create();
        vao.bind();

        glVertexAttribPointer(m_basicProgram_position, 3, GL_FLOAT, GL_FALSE, 6*sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glDrawArrays(GL_TRIANGLES, 0, len);

        glDisableVertexAttribArray(0);

        vao.release();

        vbo.release();
    }

    unsigned char res[4];
    GLint viewport[4];

    int qtx = position.x;
    int qty = position.y;

    glGetIntegerv(GL_VIEWPORT, viewport);
    int glx = viewport[2] * ((float)qtx / (float)this->width());
    int gly = viewport[3] * ((float)qty / (float)this->height());
    glReadPixels(glx, viewport[3] - gly, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, &res);

    int clickedColor = 0;
    clickedColor |= (res[0] << 16);
    clickedColor |= (res[1] << 8);
    clickedColor |= res[2];

    m_clickProgram->release();

    if (colorCodes.contains(clickedColor))
        emit regionClicked(colorCodes[clickedColor]);
}

void MeshWidget::mousePressEvent(QMouseEvent* event)
{
    this->m_clicked = true;
    this->m_mouse_pos = {event->pos().x(), event->pos().y()};
    this->update();
}

QList<QString> MeshWidget::regions()
{
    return this->m_regions.keys();
}

void MeshWidget::changeRegionState(const QString& region)
{
    bool _isActive = this->m_regions[region].isActive;
    this->m_regions[region].isActive = !_isActive;

    this->update();
}

QVector3D MeshWidget::pickColor(float angle)
{
    float theta = angle;
    float r;
    float b;
    float g;

//    while (theta < 0)
//        theta += 360;

//    while (theta >= 360)
//        theta -= 360;

    if (theta < 120) {
        g = theta / 120;
        r = 1 - g;
        b = 0;
    } else if (theta < 240) {
        b = (theta - 120) / 120;
        g = 1 - b;
        r = 0;
    } else {
        r = (theta - 240) / 120;
        b = 1 - r;
        g = 0;
    }

    return QVector3D(r, g, b);
}
