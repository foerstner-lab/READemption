test:
	pytest

coverage:
	python3-coverage run tests/test_all.py
	python3-coverage report

package:
	python3 setup.py sdist bdist_wheel
	rm -rf READemption.egg-info
	ls dist/*

clean:
	find -name "*pyc" -exec rm {} \;
	rm -rf build dist READemption.egg-info *~

build:
	 python3 setup.py bdist

package_to_pypi:
	python3 -m twine upload  --verbose  dist/*

html_doc:
	cd docs && make html && cd ..

show_html_docs:
	firefox docs/build/html/index.html &

pylint:
	pylint bin/reademption reademptionlib/* tests/*

new_release:
	@echo "* Create/checkout a release branch"
	@echo "  git branch release_v0.4.X"
	@echo "  git checkout release_v0.4.X"
	@echo "* Change bin/reademption"
	@echo "* Change setup.py"
	@echo "* Change docs/source/conf.py"
	@echo "* Change CHANGELOG.txt"
	@echo "* Create new docs"
	@echo "  make html_doc"
	@echo "* Test package creation"
	@echo "* Test doc creation"
	@echo "* make package_to_pypi"
	@echo "* git add CHANGELOG.txt bin/reademption docs/source/conf.py setup.py"
	@echo "* Commit changes e.g. 'git commit -m \"Set version to 0.4.X\"'"
	@echo "* Tag the commit e.g. 'git tag -a v0.4.X -m \"version v0.4.X\"'"
	@echo "* Merge release into dev and master"
	@echo "* Push it to github: git push"
	@echo "* Generate a new release based on this tag at"
	@echo "  https://github.com/foerstner-lab/READemption/releases/new"
	# @echo "* Upload new docs using 'make upload_doc'"

# Todo: docker_images: docker build -t konradfoerstner/reademption:0.4.1 -t konradfoerstner/reademption:latest .
