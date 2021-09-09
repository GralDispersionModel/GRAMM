# Contributing to the prognostic, non-hydrostatic mesoscale model GRAMM
Thank you very much for developing GRAMM further or for fixing bugs, so that the entire community can benefit from it!

Do not hesitate to contact the project administrators at the beginning of your work. Changes to GRAMM have to go through a complex validation
 process, so it is advantageous if the validation body, which is currently Dietmar Oettl, knows exactly the content of the changes.

## Branch Configuration

```
-- master  : production and bug fixes
-- V2XXX   : release ready commits and bug fixes
-- features/feature-xx: always branch from develop and delete after merging to develop
```

- *master* branch is inteded for production release. Keep it simple and easy to rollback
- *V2XXX*  branch is for release preparation. Only for release ready commits.


## Recommended Process

If you're developing a **new feature**

1. Create a feature branch from `V2XXX` branch
2. Branch name dependend on your new `feature`
3. When your code is ready for release, pull request to the `V2XXX` branch
4. Delete the feature branch


If you're making a **bug fix**

1. Pull request to the `V2XXX` branch
2. Add an issue tag in the commit message or pull request message

If you're making a **hot fix**, which has to be deployed immediately.
1. Pull request to `V2XXX` **and** `master` branch

## I don't want to contribute, I just have a question!
Support is provided by the [Technical University of Graz, Austria](http://lampz.tugraz.at/~gral/). 

## Found a Bug?
If you find a bug in the source code, you can help us by submitting an issue to our GitHub Repository. Even better, you can submit a Pull Request with a fix or send us an E Mail.
Please test the bug fix by one ore more projects and document the changes.

## What should I know before I get started?
GRAMM is developed on .NetCore 3.1. You can use Visual Studio or Visual Studio Code for development across platforms or Visual Studio 2019 in Windows.<br>
The released GRAMM application was compiled with GDAL and ECMWF coupling. If you want to compile without ECMWF coupling delete the constant `_ECMWF_` and remove the dependencies for GDAL in the compiler settings or in the file Source.csproj.<br>
If you want to compile with ECMWF coupling, keep the flag `_ECMWF_` and install the dependencies for GDAL. <br>
The program version including GDAL and ECMWF coupling was not tested at Linux.


## Design Decisions
For performance reasons, static jagged arrays and as few classes as possible are used (avoidance of boxing/unboxing). 

## Styleguides
We follow the Microsoft design rules.

### Git Commit Messages
* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Change array a[] to..." not "Changes array a[] to...")
* Reference issues and pull requests liberally after the first line

