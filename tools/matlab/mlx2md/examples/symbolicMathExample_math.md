# Symbolic Math Example


Copyright 2020 The MathWorks, Inc.



```matlab:Code
syms theta
f = sin(theta)
```

f = 

  
$$
 \sin \left(\theta \right)
$$

```matlab:Code
df = diff(f,theta)
```

df = 

  
$$
 \cos \left(\theta \right)
$$

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

  
$$
 \left(\begin{array}{cccc}
\cos \left(\theta \right) & -\sin \left(\theta \right) & 0 & 0\\
\sin \left(\theta \right) & \cos \left(\theta \right) & 0 & 0\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{array}\right)
$$

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

  
$$
 \left(\begin{array}{cccc}
\cos \left(\theta \right) & 0 & \sin \left(\theta \right) & 0\\
\sin \left(\theta \right) & 0 & -\cos \left(\theta \right) & 0\\
0 & 1 & 0 & 0\\
0 & 0 & 0 & 1
\end{array}\right)
$$
