
This code documentation is generated using sphinx
(http://www.sphinx-doc.org/en/master/index.html).

The following steps were taken to generate the documentation:

    mkdir rst html
    sphinx-quickstart (mostly default options, y for autodoc, y for links to
        source code, y to making .nojekyll file)
    modify conf.py file to include amoebaelib and amoebae directories in path
    change the html_theme variable to 'default'
    make html
    cd rst
    sphinx-apidoc -o ../rst/ ../../../amoebaelib/
    cp modules.rst index.rst
    cp ../conf.py .
    mkdir _static
    cd ..
    sphinx-build -b html ./rst ./html/
    open html/index.html


The following steps are required to update the documentation

    cd rst
    rm -rf * (***Careful!***)
    mkdir _static
    cp ../conf.py .
    sphinx-apidoc -o ../rst/ ../../../amoebaelib/
    cp modules.rst index.rst
    cd ..
    sphinx-build -b html ./rst ./html/
    open html/index.html

    If it doesn't work, check that the amoebae_path variable in the conf.py
    script is the actual path to the root directory of the repository.
