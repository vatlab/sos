{%- extends 'full.tpl' -%}


{%- block header -%}
{{ super() }}

{%- if nb['metadata'].get('sos',{}).get('kernels',none) is not none -%}

<style type="text/css">

table {
   padding: 0;border-collapse: collapse; }
table tr {
   border-top: 1px solid #cccccc;
   background-color: white;
   margin: 0;
   padding: 0; }
table tr:nth-child(2n) {
   background-color: #f8f8f8; }
table tr th {
   font-weight: bold;
   border: 1px solid #cccccc;
   margin: 0;
   padding: 6px 13px; }
table tr td {
   border: 1px solid #cccccc;
   margin: 0;
   padding: 6px 13px; }
table tr th :first-child, table tr td :first-child {
   margin-top: 0; }
table tr th :last-child, table tr td :last-child {
   margin-bottom: 0; }

div.input {
    display: none;
}

.input_prompt {
    display: none;
}

.output_prompt {
    display: none;
}

#nextsteps {
   color: blue;
}

</style>

<script>
function toggle_source() {
  var x = document.getElementsByClassName("input");
  if (x.length == 0) return;
  function toggle_vis(o) {
    var d = o.style.display;
    o.style.display = (d == "flex" || d == "" || d == "block" | d == "-webkit-box") ? "none": "-webkit-box";
  }

  for (i = 0; i < x.length; i++) {
    toggle_vis(x[i]);
  }
   btn = document.getElementById("1")
   if (btn.textContent.indexOf("Hide") > 0) 
      btn.textContent = "Show source code of this document";
    else
      btn.textContent = "Hide source code of this document";
}
</script>

<button id="1" type="button" onclick="toggle_source();" >Hide source code of this document</button> 

</script>

<style> 

 {% for item in nb['metadata'].get('sos',{}).get('kernels',{}) %}

{%- if item[2] -%}
.lan_{{item[0]}} {background-color: {{item[2]}} !important }

{%- else -%}
.lan_{{item[0]}} {}

{%- endif -%}

{% endfor %}

</style>


{%- endif -%}

{%- endblock header -%}


{% block codecell %}

	<div class="border-box-sizing code_cell rendered">
	{{ super() }}
	</div>

{%- endblock codecell %}

