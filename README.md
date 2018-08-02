# Dales-U

This is the basic repository for using the Dales-Urban model.

## Getting Started

### Get a copy of the repository
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

git clone instructions ...

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## How keep the submodules up-to-date

The current submodules are pre-post and dales-urban in the src directory.

First set the following git shortcuts

```
git config --global status.submoduleSummary true
git config --global pull.rebase true
git config --global alias.pula 'pull --rebase --recurse-submodules'
git config --global alias.puls 'submodule update --remote --rebase'
git config --global alias.spush 'push --recurse-submodules=on-demand'
```

and check these configurations with `git config --list`.

### Push changes in the submodules and repository

Changes in the submodule need to be added directly from within it.
So go to the submodule directory and track your changes by

```
git add changedfile
git commit -m "change message"
```

#### Push changes of the submodule only
* From within the submodule do

```
git pull
git push
```

#### Push changes of the main repository only
* From the main repository (dales-u) use

```
git push
```

#### Push changes of the main repository and submodules

* From the main repository (dales-u) use

```
git spush
```

or
* Push changes of the submodule as above. 
Then go to the main repository and add and commit the newest version of the submodule.
Push the update with

```
git push
```

### Update the repository

#### To update the submodules only:

* Go to the submodule and use 

```
git pull
```

or
* From the main repository (dales-u) use

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

#### To update the main repository only:

* From the main repository (dales-u) use

```
git pull
```

#### To update both the main repository and all submodules
* From the main repository (dales-u) use

```
git pula
```

Remember that the main repository needs to know which version of the submodule to use.
If there is a new version of the submodule it needs to be added and committed in the main repository too.
(Check it with `git status`).


## Running the tests

Explain how to run the automated tests for this system

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
