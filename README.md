# 数值计算toolbox
## **python 版本**
- 3.7
## **已有功能**
### **解线性方程组的直接法**
- 高斯-若尔当消去法
- 道立特分解法
- 追赶法
### **解线性方程组的迭代法**
- 高斯-赛德尔迭代法
- 松弛法
### **非线性方程求根**
- BM
- BMZ
- 试位法
- 简单迭代法
- 埃特金法
- 牛顿迭代法
- 重根加速法
- 弦截法
- 牛顿下山法
### **插值法**
- 拉格朗日插值法
- 牛顿插值法
## **依赖库**
- matplotlib
- numpy
- scipy
### **安装方法**
&emsp;&emsp;进入requirements.txt所在目录，执行pip install -r requirements.txt
## **项目使用说明**
### **项目安装**
一、使用git：
```cmd
git clone https://github.com/black-programmer/toolbox.git
cd toolbox
python setup.py sdist
pip install dist/toolbox-1.0.tar.gz
```
二、使用pip：
```cmd
pip install toolbox
```
### 项目导入
```python
import toolbox
...
```
具体使用方法以及包含的功能参见toolbox/doc/doc.md

测试用例位于test/test.py