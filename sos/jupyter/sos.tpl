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

.dataframe_container { max-height: 400px }
.dataframe_input {
    border: 1px solid #ddd;
    margin-bottom: 5px;
}

div.input {
    display: none;
}

.hidden_output {
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

function toggle_vis(o) {
    var d = o.style.display;
    o.style.display = (d == "flex" || d == "" || d == "block" | d == "-webkit-box") ? "none": "-webkit-box";
}

function toggle_source() {
    var btn = document.getElementById("show_1");
	var hide = true;
    if (btn.textContent.indexOf("Less") > 0) {
        btn.textContent = "Show More";
    } else {
	    hide = false;
        btn.textContent = "Show Less";
	}
    var x = document.getElementsByClassName("input");
    for (i = 0; i < x.length; i++) {
        if (hide) {
		    x[i].style.display = 'none';
		} else {
		    x[i].style.display = 'flex';
		}
    }
	var x = document.getElementsByClassName("hidden_output");
    for (i = 0; i < x.length; i++) {
        if (hide) {
		    x[i].style.display = 'none';
		} else {
		    x[i].style.display = 'block';
		}
    }
}
</script>

<button id="show_1" type="button" onclick="toggle_source();" >Show More</button> 

</script>

{%- endif -%}

{%- endblock header -%}

{%- block input -%}

	{{ super() }}

{%- endblock input -%}


{% block output %}
	{%- if cell.metadata.show_output -%}
	    {{ super() }}
    {%- else -%}
	    <div class="hidden_output">
	    {{ super() }}
		</div>
   {%- endif -%}
{% endblock output %}
