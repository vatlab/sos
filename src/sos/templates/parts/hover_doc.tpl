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

let keyword_links = {
  // cm-keyword cm-strong
  'input:' : 'https://vatlab.github.io/sos-docs/doc/user_guide/input_statement.html',
  'output:' : 'https://vatlab.github.io/sos-docs/doc/user_guide/output_statement.html',
  'depends:' : 'https://vatlab.github.io/sos-docs/doc/user_guide/depebds_statement.html',

  // cm-variable cm-sos-option
  'named_output' : 'https://vatlab.github.io/sos-docs/doc/user_guide/named_output.html',
  'output_from' : 'https://vatlab.github.io/sos-docs/doc/user_guide/output_from.html',
  'for_each': 'https://vatlab.github.io/sos-docs/doc/user_guide/for_each.html',
  'trunk_size': 'https://vatlab.github.io/sos-docs/doc/user_guide/trunk_size.html',
  'expand': 'https://vatlab.github.io/sos-docs/doc/user_guide/scripts_in_sos.html#option-expand',

  // cm-builtin cm-strong
  'run:': 'https://vatlab.github.io/sos-docs/doc/user_guide/shell_actions.html',
  'sh:': 'https://vatlab.github.io/sos-docs/doc/user_guide/shell_actions.html',
  'bash': 'https://vatlab.github.io/sos-docs/doc/user_guide/shell_actions.html',
  'R:': 'https://vatlab.github.io/sos-docs/doc/user_guide/script_actions.html',
  'Python:': 'https://vatlab.github.io/sos-docs/doc/user_guide/script_actions.html',

  // cm-meta
  '%sosrun': 'https://vatlab.github.io/sos-docs/doc/user_guide/sos_in_notebook.html#magic-sosrun',
  '%run': 'https://vatlab.github.io/sos-docs/doc/user_guide/sos_in_notebook.html#magic-run',
  '%runfile': 'https://vatlab.github.io/sos-docs/doc/user_guide/sos_in_notebook.html#magic-runfile',
  '%sossave': 'https://vatlab.github.io/sos-docs/doc/user_guide/magic_sossave.html',
  '%preview': 'https://vatlab.github.io/sos-docs/doc/user_guide/magic_preview.html',
}

function visit_sos_doc(evt) {
  window.open(keyword_links[evt.target.innerText], '_blank')
}

function add_hoverdoc(){
  let elems = ['cm-keyword cm-strong', 'cm-variable cm-sos-option', 'cm-builtin cm-strong', 'cm-meta'].map(
    cls => Array.from(document.getElementsByClassName(cls))).reduce((r, a) => r.concat(a), [])
  Array.from(elems).filter(elem => elem.innerText in keyword_links).forEach(x => {
    x.classList.add('sos_hover_doc');
    x.addEventListener('click', visit_sos_doc);
  })

}

</script>
{% endmacro %}
