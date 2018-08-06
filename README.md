# Dales-U

This is a starter kit for using the Dales-Urban model. It has the basic repository (dales-u) that contains everything you need to run Dales-U. Within it are the two submodules "pre-post" and "dales-urban" in the src directory. These libraries are currently under development and are repositories themselves. In [How keep the submodules up to date](https://gitlab.com/bss116/dales-u#how-to-keep-the-submodules-up-to-date) you find instructions on how to make sure you always have the latest version of the libraries.


## Getting Started

### Get a copy of the repository

* To get a copy of Dales-U just for your local machine go to the directory you want to store it in and use

```
git clone --recurse-submodules git@gitlab.com:dales-urbanists/dales-u.git
```

* To be able to push your changes to a remote repository, fork the project on GitLab: [Dales-U](https://gitlab.com/dales-urbanists/dales-u). 

Go to your fork, copy the ssh link (`git@gitlab.com:USERNAME/REPOSITORY.git`) and copy the project to your local machine by using

```
git clone --recurse-submodules git@gitlab.com:USERNAME/REPOSITORY.git
```

Then go to all submodules and checkout the branch you want to work on, e.g.

```
git checkout master
```

To keep your repository of Dales-U in sync with the original, check out how to [sync a fork](https://help.github.com/articles/syncing-a-fork).

### Prerequisites

Dales requires several packages installed on your local machine. The packages are gcc, gfortran, make, netcdf, open-mpi and fftw.

#### Prerequisites on high performace clusters

Load the required packages by using

```
module load packagename
```

The file "dalesmodules" in the utils directory provides a list of all required packages.

#### Prerequisites on macOS

We recommend installing the packages via [homebrew](https://brew.sh/). If you do not have homebrew, first check that XCode is installed by typing

```
xcode-select -p
``` 

If you get an error, install it with

```
xcode-select --install
```

Then install homebrew with 

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

```
To test whether homebrew is installed type

```
brew doctor
```

Install the packages netcdf (also installs gcc and gfortran), open-mpi, fftw and make (should be already installed) with

```
brew install packagename
```

You may also want to install the package nco, which is required for postprocessing in "pre-post", and ncview for easy visualisations of Dales-U output.

Homebrew installs packages in "/usr/local/" directory, make sure this path (or whereever your packages are stored) is added to the Makefile as INCDIRS and/or LIBDIRS.


### Installing

To get an executable version of Dales-U go to the directory "src/dales-urban" and compile by typing `make`.

Don't forget to run a `make clean` after making changes to the source code or Makefile and before recompiling.

## How to keep the submodules up-to-date

The current submodules are "pre-post" and "dales-urban" in the src directory.

First set the following git shortcuts

```
git config --global status.submoduleSummary true
git config --global pull.rebase true
git config --global alias.pula 'pull --rebase --recurse-submodules'
git config --global alias.puls 'submodule update --remote --rebase'
git config --global alias.spush 'push --recurse-submodules=on-demand'
```

and check these configurations with `git config --list`.

### Update the repository

#### To update the submodules only

* Go to the submodule and use 

```
git pull
```

or

* From the main repository (Dales-U) use

```
git puls
``` 

to get latest version of all submodules or just use

```
git puls modulename
``` 

to update the submodule modulename.

* To just check which updates are available you can go to the submodule directory and use

```
git fetch
git status
```

#### To update the main repository only

* From the main repository use

```
git pull
```

#### To update both the main repository and all submodules

* From the main repository use

```
git pula
```

### Publish changes in the repository

Changes to the submodule need to be reported in two steps.
1) The submodule needs to be updated.
2) The main repository (Dales-U) needs to be told to use the updated version of the submodule.

#### Push changes of the submodule

Changes to the submodule need to be added directly from within it.
Go to the submodule directory and track your changes with

```
git add changedfile
git commit -m "change message"
```
then use

```
git pull
git push
```

to update the submodule.
If you want to tell the main repository to update to this version of the submodule, go to the main repository and add the newest changes with

```
git add submodule
git commit -m "Updated submodule to newest version"
```

and then push the change with

```
git pull
git push
```

This is the recommended workflow because it is the easiest way to separate changes to the library (e.g. the Dales source code) and changes to the repository (e.g. tracking new experiment setups).

#### Push changes of the main repository only

If you make changes in the main repository and the submodule, but you do not want to update to the newest version of the submodule yet (e.g. because the changes are unstable), you can just push your changes of the main repository using

```
git pull
git push
```

from the main repository. It will then update the files but switch to the newest version of the submodule unless you commit these changes too (see above).

#### Push changes of the main repository and submodules

When you have committed changes in both the main repository and the submodules you can push them at once by using

```
git pula
git spush
```

from the main repository.
Remember that the main repository needs to know which version of the submodule to use.
If you want to update the repository to a new version of the submodule it needs to be added and committed in the repository too. (Check it with `git status`).


## Running the tests

**_Template needs to be completed from here_**

Explain how to run the automated tests for this system

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
