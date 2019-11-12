# Mackey
This is a project devoted to numerically computing the RO(G) homology of a point. 

There are three ways you can use it:

* For a quick demonstration in the case of G=C4, you can get the binary for your OS from <a href="https://github.com/NickG-Math/Mackey/tree/master/bin">here</a>.
* If you are familiar with MATLAB, then you can get an earlier prototype of this project written in MATLAB from the sister <a href="https://github.com/NickG-Math/C4-Homology">repository</a>. While it's fully documented, it should be noted that the MATLAB code is somewhat out-of-date, compared to the current project. 
* This project is a header-only library written in C++. It's more general, more modular and extensible, and about 20-30 times faster than the MATLAB version. To get started with a user-guide/tutorial, please refer to this <a href="https://nickg-math.github.io/Mackey/html/index.html">page</a>. There you will also find extensive documentation for every method and class of the project.

You can also view a graph created by this library and drawn by graphviz  <a href="https://github.com/NickG-Math/Mackey/blob/master/Multiplication_Graph.svg">here</a> (first download it and then open the svg via a browser).

What follows is a very brief installation guide taken from the more extensive <a href="https://nickg-math.github.io/Mackey/html/index.html">documentation</a>.

# Requirements
* C++17 and the standard library.
* <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a>, a header only library for matrix manipulation. I've tested this project with Eigen 3.3.7
* Optional: For improved performance you can use the Intel MKL with Eigen and further combine with OpenMP for multithreading.
* Optional: To draw the multiplication graphs you will need Graphviz.

# Installation
* To install simply clone/download this repository and include the entire folder in your path (the docs and bin subfolders are optional and not part of the source code). You will also need to do the same with Eigen.
* See this <a href="https://nickg-math.github.io/Mackey/html/use.html">page</a> for details on how to set up and call the library from your source code.
* As for compiler support, I have tested the code with the following C++ compilers: GCC 9.2 (Linux), Clang 10 (Linux and MacOS), Intel Compiler 19 (Linux and Windows), MSVC 19 (Windows). Make sure you use the option ```-std=c++17```. For more information on compiler options, see the <a href="https://nickg-math.github.io/Mackey/html/perf.html">performance</a> page.
