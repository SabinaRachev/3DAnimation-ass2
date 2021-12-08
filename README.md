ass1 :
The following files were the ones that we had modified to support the simplification algorithm:
1.Sandbox.cpp / Sandbox.h
2.inputManager.cpp
3.Viewer.cpp /Viewer.h - added a vitual function that is called initData in order to support simplification for loaded objects
4.ImGuiMenu - added the command viewer->initData() to the function that handles loading objects.

ass2:
The following files were the ones that we had modified to support the collision detection algorithm:
1.Sandbox.cpp / Sandbox.h
2.inputManager.cpp


in order to run the project:
1.clone the git repository
2.Run Cmake gui. choose the project folder and destination folder for the cpp project files.
3.Search for "texture_image"on viewerdata.cpp and replace "snake1" with the correct path of your directory.
4.Confirm correct configurations paths.
5.Copy configuration.txt from tutorial/sandBox to build/tutorial/sandBox sandBox as a startup project and compile the project (it could take few minutes);
6.Build the project
