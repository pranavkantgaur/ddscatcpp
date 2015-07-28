# DDscat.C++ #

Discrete Dipole Scattering C++ code (DDscat.C++) is a C++ version of the Fortran code DDSCAT for calculating scattering and absorption of light by irregular particles and periodic arrangement of irregular particles.

Current C++ version of the code is C++ porting of Fortran DDSCAT.

The Fortran code has been jointly developed by
<ul><li><a href='http://en.wikipedia.org/wiki/Bruce_T._Draine'>Bruce T. Draine</a> (Dept. of Astrophysical Sciences, Princeton University, draine@astro.princeton.edu) and  </li><li>Piotr J. Flatau (Scripps Institution of Oceanography, University of California San Diego, pflatau@ucsd.edu)   </li></ul>

and hosts at <a href='http://code.google.com/p/ddscat/'><pre><code>http://code.google.com/p/ddscat/</code></pre></a>

The C++ code was developed by
<ul>
<li>
Vasyl Choliy, Dr. (Kyiv Shevchenko University, Astronomy and Space Physics Department, Choliy.Vasyl(at)gmail.com), Kyiv, Ukraine.<br>
</li>
</ul>

The C++ code follows as closely as possible to the Fortran versions, starting from DDSCAT 7.3.0. Versions of C++ are marked likewise the Fortran ones, but ++ is added to the end of name string. Namely, C++ version of Fortran code 7.3.0 is DDscat.C++ 7.3.0.

Not always it is possible to map Fortran features to C++, that is why DDscat.C++ uses C++ features preferably.

DDscat.C++ 7.3.0 is a branch of Open Source project DDSCAT.

DDscat.C++ is subject to the GNU General Public License. You may copy, distribute, and/or modify the software identified as under this agreement. If you distribute copies of this software, you must give the recipients all the rights that you have.

It also means that manipulations with the DDscat.C++ are explained in Fortran version User guide. Please, refer to that documents first. New features and C++ specific issues are explained in DDscat.C++ User guide which we will refer to as 'User guide ++'.

DDscat.C++ 7.3.0 is publicly available. The DDscat.C++ 7.3.0 is then considered to be the standard version of DDscat.C++. If you choose to use it, please send email to Choliy.Vasyl(at)gmail.com  ''registering'' as a user; registered users of DDscat.C++ will be notified when updates to the code are made.

The DDscat.C++ was written as a teacher-and-student project in IT-and-physics with the intention to be useful only on the basis of Fortran code. We made a serious efforts to make the C++ code free of errors. Anyway, please, let us know if you have troubles using DDscat.C++, or if you catch some bugs in it or in User guide++.

If you use DDscat.C++ in you work, pleace, make the citation of the article

```
V.Choliy
The discrete dipole approximation code DDscat.C++: features, limitations and plans
2013, Adv.Astron.Spa.Physics, 3, 63-70.
```

# Release notes #
The code package contains:
  * h and cpp files of all components;
  * make files to build the code under Linux;
  * project files for Qt Designer with Qt 4.7,4 to build the code under Windows, Linux and Mac OS X;
  * XCode 3.0 project files to build the code under Mac OS 10.5;
  * files necessary to run examples, identical to those in Fortran versions;

The package doesn't contain any extra binary files inside. Third party xml libraries should be added by user.

The code uses the plugin paradigm. It means that the code consists of components with strictly defined communication protocol and lightweight replacement.

The code consist of:
  * FourierLib - classes and proxies to do Fourier transforms (dynamical linked library);
  * SolverLib - classes and proxies to solve equations (dynamical linked library);
  * TargetLib - classes to manage targets (dynamical linked library);
  * DDscat.C++ itself to be compiled into executable file; it consists of
    1. GreenFunctionManager
    1. DielectricManager
    1. ScatterManager
    1. TimerManager
    1. a lot of another code;
  * Readnf1 (extract field data) and Readnf2 (make a cutting of the field along the line);
  * VtrConvert;
  * Calltarget (uses TargetLib);
  * The test library.

Some features of DDscat.C++.
  * DDscat.C++ allocates the memory once, after loading the target or in the process of loading. Preliminary memory allocation is then not necessary and is omitted. Despite of that the ddscat.par left unchanged, the 8th and 9th lines of ddscat.par are read but skipped.
  * All official DDscat examples run fine. There is Tests directory with the command files in it to run all tests.
  * Readnf and VtrConverts works fine and produce the same output as Fortran version.
  * Intel MKL is not tested due to License limitings.
  * All par files have their identical xml versions. Xml parameters file is allowed, but normal text parameters file is preferable. If run without command line parameters the code analyzes if ddscat.par is present in working directory and uses it, otherwise ddscatpar.xml is searched for and used instead. The code automatically recognize what kind of parameters file is given in command line if any.
  * Nearfield produces normal text output but can produce xml output too in VTKfile format.
  * Readnf is subdivided into two programs: Readnf2 to output the field along specified line and Readnf1 to output near-target field. Readnf can use and is capable to recognize text and xml versions of readnf.par.
  * Calltarget is able to read its parameters in xml format.
  * MPI and OpenMP code are temporary switched off. See todo list below.
  * SolverLib and FourierLib are designed as plugins. They allow the usage of any external Solve or FFT engine.
  * SolverLib and FourierLib are capable to work on GPU with CUDA installed (DDscat.C++ 7.3.1 as a minimum).

# Downloading the DDscat 7.3.0++ Code and Documentation #
Access to the code tree is possible only via svn. Please, go to Source section of the site to access the code.
Users may download Delivery package from my personal page
```
http://space.univ.kiev.ua/Choliy/DDscatcpp/
```

User guide of DDscat.C++ is in Doc subfolder of main code tree or here:
```
http://space.univ.kiev.ua/Choliy/DDscatcpp/DDScatcpp.7.3.0.Userguide.pdf
```

# Compiling and building DDscat.C++ #
We made the serious efforts to make the code compilable at the famous platforms.
These are Windows (2000 and later), Linux (Ubuntu 10), Mac (OS X.10.5 and later).
Main development was done under Qt 4.7.4 at Windows XP.
The sponsor (want to stay anonymous) was so kind to allow testing of the final versions under VS2003.
Please, read DDscat.C++ User guide to start using the code.

# Acknowledgements #
Project owner wish to thanks to Bruce Draine and Piotr Flatau for their warm attitude
to DDscat.C++.