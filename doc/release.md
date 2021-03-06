# Releasing a new version of sourmash

These are adapted from the khmer release docs, originally written by
Michael Crusoe.

Remember to update release numbers/RC in:

* this document

## Testing a release

0\. First things first: check if Read the Docs is building properly for master.
The build on [Read the Docs] should be passing,
and also check if the [rendered docs] are updated.

[Read the Docs]: https://readthedocs.org/projects/sourmash/builds/
[rendered docs]: https://sourmash.readthedocs.io/en/latest/

1\. The below should be done in a clean checkout:
```
cd $(mktemp -d)
git clone git@github.com:dib-lab/sourmash.git
cd sourmash
```

2\. Set your new version number and release candidate
(you might want to check [the releases page] for next version number):
```
new_version=3.0.0
rc=rc1
```
and then tag the release candidate with the new version number prefixed by the letter 'v':
```
git tag -a v${new_version}${rc}
git push --tags git@github.com:dib-lab/sourmash.git
```
[the releases page]: https://github.com/dib-lab/sourmash/releases

3\. Test the release candidate. Bonus: repeat on macOS:
```
set -e
python -m pip install -U setuptools pip virtualenv wheel

cd ..
python -m venv testenv1
python -m venv testenv2
python -m venv testenv3
python -m venv testenv4

# First we test the tag

cd testenv1
source bin/activate
git clone --depth 1 --branch v${new_version}${rc} https://github.com/dib-lab/sourmash.git
cd sourmash
python -m pip install -U setuptools pip
python -m pip install -r requirements.txt
make test

# Secondly we test via pip

cd ../../testenv2
deactivate
source bin/activate
python -m pip install -U setuptools pip
python -m pip install -e git+https://github.com/dib-lab/sourmash.git@v${new_version}${rc}#egg=sourmash[test]
cd src/sourmash
make test
make dist
cp dist/sourmash*tar.gz ../../../testenv3/

# Is the distribution in testenv2 complete enough to build another
# functional distribution?

cd ../../../testenv3/
deactivate
source bin/activate
python -m pip install -U setuptools pip
python -m pip install sourmash*tar.gz
python -m pip install pytest
tar xzf sourmash-${new_version}${rc}.tar.gz
cd sourmash-${new_version}${rc}
python -m pip install -r requirements.txt
make dist
cp -a ../../sourmash/tests/test-data tests/  ## We don't ship the test data, so let's copy it here
make test
```

4\. Publish the new release on the testing PyPI server.
You will need to [change your PyPI credentials].
We will be using `twine` to upload the package to TestPyPI and verify
everything works before sending it to PyPI:
```
python -m pip install twine
twine upload --repository-url https://test.pypi.org/legacy/ dist/sourmash-${new_version}${rc}.tar.gz
```
Test the PyPI release in a new virtualenv:
```
cd ../../testenv4
deactivate
source bin/activate
python -m pip install -U setuptools pip wheel
# install as much as possible from non-test server!
python -m pip install screed pytest numpy matplotlib scipy khmer ijson bam2fasta deprecation cffi
python -m pip install -i https://test.pypi.org/simple --pre sourmash
sourmash info  # should print "sourmash version ${new_version}${rc}"
```

[change your PyPI credentials]: https://packaging.python.org/tutorials/packaging-projects/#uploading-the-distribution-archives

5\. Do any final testing:

 * check that the binder demo notebook is up to date
 * check wheels from GitHub releases

## How to make a final release

When you've got a thoroughly tested release candidate,
cut a release like so:

1\. Create the final tag. Write the changes from previous version in the tag commit message. `git log --oneline` can be useful here, because it can be used to compare the two versions (and hopefully we used descriptive PR names and commit messages). An example comparing `2.2.0` to `2.1.0`:
`git log --oneline v2.1.0..v2.2.0`

```
cd ../sourmash
git tag -a v${new_version}
```

2\. Publish the new release on PyPI (requires an authorized account).
```
make dist
twine upload dist/sourmash-${new_version}.tar.gz
```

3\. Delete the release candidate tag and push the tag updates to GitHub:
```
git tag -d v${new_version}${rc}
git push --tags git@github.com:dib-lab/sourmash.git
git push --delete git@github.com:dib-lab/sourmash.git v${new_version}${rc}
```
4\. Add the release on GitHub, using the tag you just pushed.
Name it 'version X.Y.Z', and copy and paste in the release notes.

5\. Upload wheels from GitHub Releases to PyPI.
[`hub`](https://hub.github.com/) makes this easier,
but you can also manually download all the files from [the releases page].
```
mkdir -p wheel && cd wheel
hub release download v${new_version}
twine upload *.whl
```

## Bioconda

The BiocondaBot has an `autobump` feature that should pick up new releases from PyPI, and open a PR in Bioconda. Review any changes
(especially dependency versions, since these don't get picked up).

An example PR for [`2.1.0`](https://github.com/bioconda/bioconda-recipes/pull/17113).

## Announce it!

If a bioinformatics software is released and no one tweets, is it really released?

Examples:

- [3.1.0](https://twitter.com/luizirber/status/1217639572202409984)
- [3.0.0](https://twitter.com/luizirber/status/1213588144458649600)
- [2.3.0](https://twitter.com/luizirber/status/1198027116396171264)
- [2.2.0](https://twitter.com/luizirber/status/1179126660911661057)
- [2.1.0](https://twitter.com/luizirber/status/1166910335120314369)
- [2.0.1](https://twitter.com/luizirber/status/1136786447518711808)
- [2.0.0](https://twitter.com/luizirber/status/1108846466502520832)

## To test on a blank Ubuntu system

```
apt-cache update && apt-get -y install python-dev libfreetype6-dev && \
python -m pip install sourmash[test]
```
