#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "resource_manager.h"
#include "canvas.h"

#include <iostream>

// GLFW function declarations
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);

// dimensions of the screen
const unsigned int WINDOW_WIDTH = 800;
const unsigned int WINDOW_HEIGHT = 600;

Canvas* canvas;

void APIENTRY glDebugOutput(GLenum source, 
    GLenum type, 
    unsigned int id, 
    GLenum severity, 
    GLsizei length, 
    const char *message, 
    const void *userParam)
{
    // ignore non-significant error/warning codes
    if(id == 131169 || id == 131185 || id == 131218 || id == 131204) return; 

       std::cout << "---------------" << std::endl;
       std::cout << "OpenGL debug message (" << id << "): " <<  message << std::endl;
}


int main () {
    // create window
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, true);
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif
    glfwWindowHint(GLFW_RESIZABLE, false);

    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Nomenclatura", nullptr, nullptr);
    glfwMakeContextCurrent(window);

    // glad: load OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glfwSetKeyCallback(window, key_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // OpenGL configuration
    glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS); 
    glDebugMessageCallback(glDebugOutput, nullptr);
    glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr, GL_TRUE);

    unsigned int cursorType = GLFW_CROSSHAIR_CURSOR;
    GLFWcursor* cursor = glfwCreateStandardCursor(cursorType);
    glfwSetCursor(window, cursor);

    canvas = new Canvas(glm::vec2(WINDOW_WIDTH, WINDOW_HEIGHT));
    canvas->init();

    // deltaTime variables
    float deltaTime = 0.0f;
    float lastFrame = 0.0f;
    
    while (!glfwWindowShouldClose(window)) {
        // calculate deltaTime
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        glfwPollEvents();

        // manage user input
        for (int i = 0; i < 8; i++) {
            if (glfwGetMouseButton(window, i) == GLFW_PRESS) {
                canvas->mouseKeys[i] = true;
            }
            else if (glfwGetMouseButton(window, i) == GLFW_RELEASE) {
                canvas->mouseKeys[i] = false;
            }
        }
        static double mousePosX, mousePosY;
        glfwGetCursorPos(window, &mousePosX, &mousePosY);
        canvas->mousePos = glm::vec2(mousePosX, mousePosY);
        canvas->processInput();

        // render
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        if (canvas->cursor != cursorType) {
            cursorType = canvas->cursor;
            cursor = glfwCreateStandardCursor(cursorType);
            glfwSetCursor(window, cursor);
        }

        canvas->render();

        glfwSwapBuffers(window);
    }
    
    ResourceManager::Clear();
    glfwTerminate();



    return 0;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (key >= 0 && key < 1024) {
        if (action == GLFW_PRESS)
            canvas->keys[key] = true;

        else if (action == GLFW_RELEASE)
            canvas->keys[key] = false;
    }
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
}