# README

This is a fork of the excellent [manuscript template from the Poisot Lab](https://github.com/PoisotLab/manuscript-template), which takes a manuscript 
written in Markdown (`manuscript.md`) and builds an HTML version hosted on Github pages and produces preprint and copyediting versions of the manuscript as pdfs, rendered with TeX. This fork has all my css stuff, etc. 

## How to use the manuscript template?

- Click on the "Use this template" button
- Go to your Github profile > `Settings` > `Developer Settings` > `Personal access tokens` and create and copy a new token. Name it after the generated repo, and give it 
- Go to your generated repositories settings > `Secrets` > `New secret`. Call the new secret `ACCESS_TOKEN` and paste the value from the token created in the previous step  
- Trigger a build by pushing a commit
- Go to `http://you.github.io/repo-name/` to view the html version, PDF, and LaTeX source
- Edit this README!
- Add your references to the `references.bib` file
- Edit the `metadata.json` file to add the title, abstract, authors
- Edit the `manuscript.md` file to make changes to your text

## What are things to *not* do?

- Do not add or remove fields in `metadata.json`
- Do not modify the version numbers in the build matrix of the github workflow

## Can I use this to build the pdf locally?

- Yes, but you're on your own
