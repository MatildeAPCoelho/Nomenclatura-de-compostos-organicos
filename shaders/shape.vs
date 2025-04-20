#version 330 core
layout (location = 0) in vec2 vertex;

uniform mat4 model;

uniform vec2 windowSize = vec2(1, 1);

void main() {
    vec4 pos = model * vec4(vertex.xy, 0.0, 1.0);
    // convert from window coordinates
    pos.x = pos.x * 2.0f/windowSize.x - 1.0f;
    pos.y = pos.y * -2.0f/windowSize.y + 1.0f;
    

    gl_Position = pos;
}