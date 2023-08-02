CXX=g++
CXXFLAGS=-IImGui -Wall -std=c++14 -fexceptions

OBJECTS=main.o\
	ImGui/imgui.o\
	ImGui/imgui_draw.o\
	ImGui/imgui-SFML.o\
	ImGui/imgui_tables.o\
	ImGui/imgui_widgets.o

run: ${OBJECTS}
	${CXX} -o run -lGL -lsfml-graphics -lsfml-audio -lsfml-network -lsfml-system -lsfml-window ${OBJECTS}

clean:
	rm -f ${OBJECTS}


