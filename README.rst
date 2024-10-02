POLSYS-RS
=========

This crate provides a Rust wrapper around the ``POLSYS_PLP`` Fortran library
by Wise, Sommese, and Watson in their `ACM paper <https://doi.org/10.1145/347837.347885>`__.
The library can solve a system of N polynomial equations with complex
coefficients with N unknowns using a globally convergence homotopy method.

This is currently **very experimental** since doing Fortran bindings is hard!

Example
=======

TODO

Development
===========

Links

* `Code <https://github.com/alexfikl/polsys-rs>`__.

Other known libraries that may do a better job if they had Rust bindings:

* `pypolsys <https://github.com/nennigb/pypolsys>`__ (Python). This has been a
  major inspiration for this wrapper, since we have similar goals.
* `PHCpack <https://github.com/janverschelde/PHCpack>`__ (Ada, C, Julia, Maple, Python).
* `Bertini 2.0 <https://github.com/bertiniteam/b2>`__ (C++, Python).
* `HomotopyContinuation.jl <https://www.juliahomotopycontinuation.org/>`__ (Julia).

License
=======

The Rust and Fortran code in this library is MIT licensed. However, the Fortran code
for the ``POLSYS_PLP`` library is under the `ACM Software License Agreement
<https://www.acm.org/publications/policies/software-copyright-notice>`__. If you
want to use this library for anything serious, make sure to enquire with the
original authors.
