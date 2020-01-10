import setuptools

with open("README.md", "r", encoding="gbk") as fh:
    long_description = fh.read()

setuptools.setup(
    name="toolbox", # Replace with your own username
    version="1.0",
    author="black-programmer",
    author_email="2860330441@qq.com",
    description=u"数值计算toolbox",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/black-programmer/toolbox",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
