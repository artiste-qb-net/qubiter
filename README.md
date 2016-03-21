#Qubiter at GitHub#
##What is Qubiter?##

The Qubiter project aims to eventually provide a full suite of tools, written mostly in Python, for designing and simulating quantum circuits on classical computers. (So it will address the needs of gate model, not annealer, quantum computer engineers only)

An earlier C++ computer program also called Qubiter (see http://www.ar-tiste.com/qubiter.html), written by Robert R. Tucci, did only quantum compiling. This newer project will eventually include a quantum compiler similar to earlier Qubiter, but one written in Python. But this new project will also include much more than that.

On this initial release, we have included classes for reading and writing quantum circuit files. Also for embedding a circuit inside a larger one. And last but not least, we've included a simulator.

The simulator hasn't been bench-marked but should be pretty fast, since it relies on Numpy, which is a Python wrapper for C code.

We use the quantum Fourier Transform circuit as an example of how to write quantum circuit files. Plus each class has a main method at the end giving more examples.

The quantum circuits are stored as text files, which allows easy exchange between QC engineers.

The quantum circuits are draw in ASCII (not in postscript or a proprietary format). We hope to convince you that ASCII drawings of quantum circuits are surprisingly beautiful, expressive, and convenient, really all you need, plus, unlike other formats, they are super easy to edit. Using other formats might require you to master something difficult like postscript in order to write/edit circuit files. This is totally unnecessary!

Quantum Fog at GitHub is a twin project started by the same people. We hope that eventually Quantum Fog will call Qubiter to perform some tasks, like quantum compiling and simulating.

Qubiter at GitHub is licensed under the BSD license (3 clause version) with an added clause at the end, taken almost verbatim from the Apache 2.0 license, granting additional Patent rights. See `LICENSE.md`.

##Contributors##

(Alphabetical Order)
* Dekant, Henning
* Tucci, Robert


