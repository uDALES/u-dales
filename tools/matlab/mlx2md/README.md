# Converting MATLAB Live Script to Markdown
Copyright 2020 The MathWorks, Inc.
[![View livescript2markdown: MATLAB's live scripts to markdown.  on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://jp.mathworks.com/matlabcentral/fileexchange/73993-livescript2markdown-matlab-s-live-scripts-to-markdown)

## Imporant Note (2023/9/14)

As of R2023b, you can use [export](https://jp.mathworks.com/help/matlab/ref/export.html) function of MATLAB to export markdown from livescript.


## NOTE (2020/02/10)

When exporting to LaTeX right after running the livescript, it's observed that the figures will be exported as eps files or not at all
if the livescript contains more than 20 figures.

I suggest that **you close the script and reopen and then export to latex.**

# Introduction

This repository provides a functions to convert your live scripts to markdown file. I hope this function makes your life easy to document your repository.

[English instruction](doc/README_EN.md)

[日本語はこちら](doc/README_JP.md)

I've checked the function with multiple live scripts but please note that it's not perfect. It's expected that you need some manual editing.



  
***
### Feedback

Hope it accelerates your MATLAB life. Any comment and suggestions are always appreciated.


