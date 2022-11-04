## Description

Please include a summary of the change and which issue is fixed.  
Please also include relevant motivation and context.  
List any dependencies that are required for this change.

Fixes # (issue)

## Type of change

Please delete options that are not relevant.

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation updated

## Conventional Commits guidelines

- [ ] I made sure the PR title follows the 
https://www.conventionalcommits.org/en/v1.0.0/


Changes to workflow inputs (sample table and/or configs)
* major (add **BREAKING CHANGE:** in the beginning of the PR title)
    * more fields/properties are required 
    * existing ones are dropped entirely 
* minor (add **feat:** in the beginning of the PR title)
    * optional fields/properties are added
    * required ones are made optional

Changes to workflow outputs
* major (add **BREAKING CHANGE:** in the beginning of the PR title)
    * changes lead to removal/change of existing outputs (format or location)
* minor (add **feat:** in the beginning of the PR title)
    * additional outputs are generated
    * content (but not format or location) of existing outputs changes

Everything else: patch
(add any other conventional commit in the beginning of the PR title)


## Checklist:

- [ ] My code changes follow the style of this project
- [ ] I have performed a self-review of my own code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing tests pass locally with my changes
- [ ] I have updated the project's documentation
