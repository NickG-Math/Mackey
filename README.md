# Mackey
This is a project devoted to numerically computing RO(G) homology. 

There are two ways you can use it:

* For a quick demonstration in the case of G=C4, you can get the binary for your OS from <a href="https://github.com/NickG-Math/Mackey/releases">here</a>. Note: The Mac and Linux binaries must be executed from a terminal.
* To get started with using the C++ library, there is a user-guide/tutorial on this <a href="https://nickg-math.github.io/Mackey/html/index.html">page</a>. There, you will also find extensive documentation for every method and class of the project.

You can also view a graph created by this library and drawn by graphviz  <a href="https://github.com/NickG-Math/Mackey/blob/master/Multiplication_Graph.svg">here</a> (first download it and then open the svg via a browser).

What follows is a very brief installation guide taken from the more extensive <a href="https://nickg-math.github.io/Mackey/html/index.html">documentation</a>.

# Requirements
* A C++17 compiler.
* <a href="http://eigen.tuxfamily.org/index.php?title=Main_Page">Eigen</a>, a header only library for matrix manipulation. This project has been tested with Eigen 3.3.9
* Optional: To serialize the results of computations, you can use the  <a href="https://uscilab.github.io/cereal/">cereal</a> library.
* Optional: To draw the multiplication graphs you will need Graphviz.

# Installation
* To install simply clone/download this repository and include the "source" folder in your path. You will also need to do the same with Eigen (and optionally, cereal).
* See this <a href="https://nickg-math.github.io/Mackey/html/use.html">page</a> for details on how to set up and call the library from your source code.
