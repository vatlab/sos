{% macro css() %}

<style type="text/css">

.sos_hover_doc:hover {
  cursor: pointer;
  background: lightgray;
  opacity: 50%;
}
</style>

{% endmacro %}

{% macro js() %}
<script>

let sos_keywords = new Set([
  // keywords
  'input:', 'output:', 'depends:', 'task:',
  // options
  'expand', 'provides',
  // input options
  'for_each', 'paired_with', 'group_with', 'pattern',
  'named_output', 'output_from', 'concurrent',
  // functions
  'sos_run', 'fail_if', 'done_if', 'warn_if', 'skip_if',
  'get_output', 'expand_pattern',
  // task options
  'walltime', 'cores', 'mem', 'queue', 'to_host', 'nodes',
  'from_host', 'trunk_size', 'trunk_workers', 'workdir',
  'shared', 'env', 'prepend_path', 'tags',
  // targets
  'file_target', 'executable', 'sos_variable', 'env_variable',
  'sos_step', 'dynamic', 'remote', 'system_resource',
  'Py_Module', 'R_Library',
  // actions
  'run:', 'script:', 'report:', 'bash:', 'sh:',
  'csh:', 'tcsh', 'zsh:', 'perl:', 'ruby:', 'node:',
  'pandoc:', 'docker_build:', 'download:', 'julia:',
  'R:', 'matlab:', 'python:', 'python2:', 'python3:',
  'singularity_build:',
  // action options
  'args', 'allow_error', 'container', 'engine', 'env',
  'stdout', 'stderr', 'tracked',
  // magics
  '%capture', '%cd', '%env',
  '%expand', '%get', '%matplotlib', '%preview',
  '%put', '%render', '%run', '%runfile', '%revisions',
  '%save', '%set', '%sessioninfo', '%shutdown', '%sosrun',
  '%convert', '%task', '%use', '%with'
  ]
)

function visit_sos_doc(evt) {
  window.open(`https://vatlab.github.io/sos-docs/redirect_doc.html?${evt.target.innerText}`, '_blank');
}

function add_hoverdoc(){
  let elems = ['cm-keyword cm-strong', 'cm-variable cm-sos-option', 'cm-builtin cm-strong', 'cm-meta'].map(
    cls => Array.from(document.getElementsByClassName(cls))).reduce((r, a) => r.concat(a), [])
  Array.from(elems).filter(elem => sos_keywords.has(elem.innerText)).forEach(x => {
    x.classList.add('sos_hover_doc');
    x.addEventListener('click', visit_sos_doc);
  })

}

</script>
{% endmacro %}
