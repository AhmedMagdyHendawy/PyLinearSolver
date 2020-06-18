import setuptools
import os
with open("README.md", "r") as fh:
    long_description = fh.read()
    
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setuptools.setup(
     name='PyLinearSolver',  
     version='0.1.0',
     author="Ahmed Hendawy, Bernhard Rolle",
     author_email="ahmedmagdyahmed1996@outlook.com",
     description="A python interface package for many well known linear solvers libraries in many languages, Julia, C++, etc...",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/AhmedMagdyHendawy/PyLinearSolver",
     packages=setuptools.find_packages(),
     install_requires=requirements,
     project_urls={
        # "Documentation": "https://docs.example.com/HelloWorld/",
        "Source Code": "https://github.com/AhmedMagdyHendawy/PyLinearSolver",
    },
    data_files = [("", ["LICENSE"])]
 )
os.system("wget https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz && tar -xvzf julia-1.4.2-linux-x86_64.tar.gz && cp -r julia-1.4.2 /opt/ && ln -s /opt/julia-1.4.2/bin/julia /usr/local/bin/julia")
os.system("python -c 'import julia; julia.install()'")