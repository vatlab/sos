{%- extends 'full.tpl' -%}


{%- block header -%}
{{ super() }}

{%- if nb['metadata'].get('sos',{}).get('kernels',none) is not none -%}

<style type="text/css">

table {
   padding: 0;
   border-collapse: collapse; }
thead {
    border-bottom-width: 1px;
    border-bottom-color: rgb(0,0,0);
    border-bottom-style: solid;
}
table tr {
   border: none;
   background-color: white;
   margin: 0;
   padding: 0; }
table tr:nth-child(2n) {
   background-color: #f8f8f8; }
table tr th {
   font-weight: bold;
   border: none;
   margin: 0;
   padding: 6px 13px; }
table tr td {
   border: none;
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

.rendered_html table {
  border: none;
}

.sos_hint {
  color: rgba(0,0,0,.4);
  font-family: monospace;
  display: none;
}

.output_stderr {
    display: none;
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

.output_area .prompt {
    display: none;
}

.output_prompt {
    display: none;
}

#nextsteps {
   color: blue;
}

.scatterplot_by_rowname div.xAxis div.tickLabel {
    transform: translateY(15px) translateX(15px) rotate(45deg);
    -ms-transform: translateY(15px) translateX(15px) rotate(45deg);
    -moz-transform: translateY(15px) translateX(15px) rotate(45deg);
    -webkit-transform: translateY(15px) translateX(15px) rotate(45deg);
    -o-transform: translateY(15px) translateX(15px) rotate(45deg);
    /*rotation-point:50% 50%;*/
    /*rotation:270deg;*/
}

div.cell {
    padding: 0pt;
    border-width: 0pt;
}
.sos_dataframe td, .sos_dataframe th, .sos_dataframe tr {
    white-space: nowrap;
    border: none;
}

.sos_dataframe tr:hover {
    background-color: #e6f2ff;
}

.display_control_panel  {
    padding: 10pt;
    left: 5px;
    top: 5px;
    position: fixed;
    z-index: 1000;
}

.display_control_panel:hover {
    background: rgb(224, 234, 241);
}
.display_checkboxes {
    margin-top: 5pt;
}
.display_control_panel:hover .display_control {
    display: block;
    opacity: 100;
}
.display_control_panel .display_control {
    opacity: 0;
}




{%- if nb['metadata'].get('sos',{}).get('kernels',none) is not none -%}

{% for item in nb['metadata'].get('sos',{}).get('kernels',{}) %}

{%- if item[2] -%}
.lan_{{item[0]}} .input_prompt { background-color: {{item[3]}} !important }  

{%- else -%}
.lan_{{item[0]}} {}

{%- endif -%}

{% endfor %}

{%- endif -%}
</style>




<script>


function filterDataFrame(id) {
    var input = document.getElementById("search_" + id);
    var filter = input.value.toUpperCase();
    var table = document.getElementById("dataframe_" + id);
    var tr = table.getElementsByTagName("tr");

    // Loop through all table rows, and hide those who don't match the search query
    for (var i = 1; i < tr.length; i++) {
        for (var j = 0; j < tr[i].cells.length; ++j) {
            var matched = false;
            if (tr[i].cells[j].innerHTML.toUpperCase().indexOf(filter) != -1) {
                tr[i].style.display = "";
                matched = true
                break;
            }
            if (!matched)
                tr[i].style.display = "none";
        } 
    }
}

function sortDataFrame(id, n, dtype) {
    var table = document.getElementById("dataframe_" + id);

    var tb = table.tBodies[0]; // use `<tbody>` to ignore `<thead>` and `<tfoot>` rows
    var tr = Array.prototype.slice.call(tb.rows, 0); // put rows into array

    if (dtype === 'numeric') {
        var fn = function(a, b) { 
            return parseFloat(a.cells[n].textContent) <= parseFloat(b.cells[n].textContent) ? -1 : 1;
        }
    } else {
        var fn = function(a, b) {
            var c = a.cells[n].textContent.trim().localeCompare(b.cells[n].textContent.trim()); 
            return c > 0 ? 1 : (c < 0 ? -1 : 0) }
    }
    var isSorted = function(array, fn) {
        if (array.length < 2)
            return 1;
        var direction = fn(array[0], array[1]); 
        for (var i = 1; i < array.length - 1; ++i) {
            var d = fn(array[i], array[i+1]);
            if (d == 0)
                continue;
            else if (direction == 0)
                direction = d;
            else if (direction != d)
                return 0;
            }
        return direction;
    }

    var sorted = isSorted(tr, fn);

    if (sorted == 1 || sorted == -1) {
        // if sorted already, reverse it
        for(var i = tr.length - 1; i >= 0; --i)
            tb.appendChild(tr[i]); // append each row in order
    } else {
        tr = tr.sort(fn);
        for(var i = 0; i < tr.length; ++i)
            tb.appendChild(tr[i]); // append each row in order
    }
}




function toggle_source() {
    var btn = document.getElementById("show_cells");
    if (btn.checked) {
        $('div.input').css('display', 'flex');
        $('.hidden_output').show();
        // this somehow does not work.
        $('div.cell').css('padding', '5pt').css('border-width', '1pt');
    } else {
        $('div.input').hide();
        $('.hidden_output').hide();
        $('div.cell').css('padding', '0pt').css('border-width', '0pt');
    }
}

function toggle_prompt() {
    var btn = document.getElementById("show_prompt");
    if (btn.checked) {
        $('.output_prompt').show();
        $('.input_prompt').show();
        $('.output_area .prompt').show();
    } else {
        $('.output_prompt').hide();
        $('.input_prompt').hide();
        $('.output_area .prompt').hide();
    }
}

function toggle_messages() {
    var btn = document.getElementById("show_messages");
    if (btn.checked) {
        $('.sos_hint').show();
        $('.output_stderr').show();
    } else {
        $('.output_stderr').hide();
        $('.sos_hint').hide();
    }
}

</script>


<div class='display_control_panel'>
    <div class="display_control">
Display content:<br>
<div class="display_checkboxes">
<input type="checkbox" id="show_cells" name="show_cells" onclick="toggle_source()">
<label for="show_cells">All cells</label>
<br>
<input type="checkbox" id="show_prompt" name="show_prompt" onclick="toggle_prompt()">
<label for="show_prompt">Prompt</label>
<br>
<input type="checkbox" id="show_messages" name="show_messages" onclick="toggle_messages()">
<label for="show_messages">Messages</label>
</div>
   </div>

</div>


</script>

{%- endif -%}

{%- endblock header -%}

{%- block input -%}

    {%- if 'scratch' in cell.metadata.tags -%}
    {%- else -%}
        {{ super() }}
   {%- endif -%}
{%- endblock input -%}


{% block output %}
    {%- if 'report_output' in cell.metadata.tags -%}
        {{ super() }}
    {%- elif 'scratch' in cell.metadata.tags -%}
    {%- else -%}
        <div class="hidden_output">
        {{ super() }}
        </div>
   {%- endif -%}
{% endblock output %}

{% block markdowncell %}
    {%- if 'hide_output' in cell.metadata.tags -%}
		<div class="hidden_output">
        {{ super() }}
		</div>
    {%- elif 'scratch' in cell.metadata.tags -%}
    {%- else -%}
        {{ super() }}
   {%- endif -%}
{%- endblock markdowncell -%}


{% block codecell %}

{%- if cell['metadata'].get('kernel',none) is not none -%}
    <div class="rendered lan_{{cell['metadata'].get('kernel', none)}}">
    {{ super() }}
    </div>
{%- else -%}
    {{ super() }}
{%- endif -%}

{%- endblock codecell %}
