<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Proposed actions](#proposed-actions)
  - [Execution of scripts](#execution-of-scripts)
  - [Utility actions](#utility-actions)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Proposed actions

### Execution of scripts

* ``run()``
* ``bash()``
* ``R()``
* ``perl()``
* ``python()``

### Utility actions

* ``fail_if()`` Terminate the execution of the pipeline if some condition is not met. This is SoS's version of  CWL's `requirements` feature.

* ``warn_if()`` Give a warning message if some condition is not met.