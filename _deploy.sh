#!/bin/sh

set -e

#[ -z "${GITHUB_PAT}" ] && exit 0
#[ "${TRAVIS_BRANCH}" != "main" ] && exit 0

#git config --global user.email "zouhua1@outlook.com"
#git config --global user.name "Hua Zou"

#git clone -b gh-pages https://${GITHUB_PAT}@github.com/${TRAVIS_REPO_SLUG}.git book-output
#cd book-output
#cp -r ../_book/* ./
git add -A
git commit -m "Update the code" || true
#git push -q origin gh-pages
git push origin main

