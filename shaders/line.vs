#version 330 core
layout (location = 0) in vec2 pos;

uniform vec2 pos1;
uniform vec2 pos2;

uniform vec2 windowSize = vec2(1, 1);

void main() {
    // line interpolation
    float xpos = pos1.x + pos.x * (pos2.x-pos1.x);
    float ypos = pos1.y + pos.y * (pos2.y-pos1.y);
    // convert from window coordinates
    xpos = xpos * 2.0f/windowSize.x - 1.0f;
    ypos = ypos * -2.0f/windowSize.y + 1.0f;
    
    gl_Position = vec4(xpos, ypos, 0.0f, 1.0f);
}