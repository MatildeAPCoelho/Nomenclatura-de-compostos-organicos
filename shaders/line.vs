#version 330 core
layout (location = 0) in vec2 pos;

uniform vec2 pos1;
uniform vec2 pos2;

uniform vec2 windowSize = vec2(1, 1);
uniform float zoom = 1.0f;
uniform vec2 offset = vec2(0, 0);

void main() {
    // line interpolation
    float xpos = pos1.x + pos.x * (pos2.x-pos1.x);
    float ypos = pos1.y + pos.y * (pos2.y-pos1.y);
    // convert from window coordinates
    xpos = ((xpos + offset.x) * 2.0f/windowSize.x - 1.0f) * zoom;
    ypos = ((ypos + offset.y) * -2.0f/windowSize.y + 1.0f) * zoom;
    
    gl_Position = vec4(xpos, ypos, 0.0f, 1.0f);
}