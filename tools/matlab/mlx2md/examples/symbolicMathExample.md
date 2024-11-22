# Symbolic Math Example


Copyright 2020 The MathWorks, Inc.



```matlab:Code
syms theta
f = sin(theta)
```

f = 

   <img src="https://latex.codecogs.com/gif.latex?&space;\sin&space;\left(\theta&space;\right)"/>

```matlab:Code
df = diff(f,theta)
```

df = 

   <img src="https://latex.codecogs.com/gif.latex?&space;\cos&space;\left(\theta&space;\right)"/>

```matlab:Code
fplot(f)
hold on
fplot(df)
hold off
```


![figure_0.png](symbolicMathExample_images/figure_0.png)


```matlab:Code
a = [cos(theta), -sin(theta),0,0
    sin(theta),cos(theta),0,0
    0,0,1,0
    0,0,0,1]
```

a = 

   <img src="https://latex.codecogs.com/gif.latex?&space;\left(\begin{array}{cccc}&space;\cos&space;\left(\theta&space;\right)&space;&&space;-\sin&space;\left(\theta&space;\right)&space;&&space;0&space;&&space;0\\&space;\sin&space;\left(\theta&space;\right)&space;&&space;\cos&space;\left(\theta&space;\right)&space;&&space;0&space;&&space;0\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\end{array}\right)"/>

```matlab:Code
b = [1,0,0,0
    0,0,-1,0
    0,1,0,0
    0,0,0,1]
```


```text:Output
b = 4x4    
     1     0     0     0
     0     0    -1     0
     0     1     0     0
     0     0     0     1

```


```matlab:Code
a*b
```

ans = 

   <img src="https://latex.codecogs.com/gif.latex?&space;\left(\begin{array}{cccc}&space;\cos&space;\left(\theta&space;\right)&space;&&space;0&space;&&space;\sin&space;\left(\theta&space;\right)&space;&&space;0\\&space;\sin&space;\left(\theta&space;\right)&space;&&space;0&space;&&space;-\cos&space;\left(\theta&space;\right)&space;&&space;0\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\end{array}\right)"/>
