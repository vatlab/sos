{% import 'parts/hover_doc.tpl' as doc %}

<!DOCTYPE html>
<html lang="en">
   <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
      <link type="text/css" rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/jquery-dropdown/2.0.3/jquery.dropdown.css" />
      <link type="text/css" rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/codemirror.css">
      <link type="text/css" rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/theme/{{theme}}.css">
      <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/css/bootstrap.min.css" integrity="sha384-WskhaSGFgHYWDcbwN70/dfYBj47jz9qbsMId/iRN3ewGhXQFZCSftd1LZCfmhktB" crossorigin="anonymous">
      <link href="https://use.fontawesome.com/releases/v5.0.6/css/all.css" rel="stylesheet" type="text/css">
      <title>{{workflow_name}}</title>
      <style type="text/css">
         {% include "parts/sos_report.css" %}
      </style>
      {{ doc.css() }}
   </head>
   <body>
      <h2 class='mt-0'>{{basename}}
         {% if url %}
         <a href="{{url}}"><i class="fas fa-external-link-alt"></i></a>
         {% endif %}
      </h2>
      <div class="file">
         <div class="fileheader">
            <div class="fileinfo">
               {{ script.splitlines() | length }} lines
               <span class="file-info-divider"></span>
               {{ script | length | filesizeformat}}
            </div>
         </div>
         <div class="filecontent">
            <textarea rows="{{ script.splitlines() | count }}" id="source-code" class="sos-source" name="code">{{ script }}</textarea>
         </div>
      </div>
      <footer>
         <div class="float-left"><a class="sos-logo" href="https://vatlab.github.io/sos-docs">
         <img src="http://vatlab.github.io/sos-docs/img/sos_icon.svg" alt="sos_icon">
         </a> &nbsp;See the output of command <samp>sos run {{basename}} -h</samp> for usage information</div>
         <div class="float-right"><a href="https://vatlab.github.io/sos-docs/">SoS</a> version <samp>{{sos_version}}</samp></div>
      </footer>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery-dropdown/2.0.3/jquery.dropdown.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js" integrity="sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49" crossorigin="anonymous"></script>
      <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/js/bootstrap.min.js" integrity="sha384-smHYKdLADwkXOn1EmN1qk/HfnUcbVRZyYmZ4qpPea6sjB/pTJ0euyQp0Mk8ck+5T" crossorigin="anonymous"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/codemirror.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/python/python.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/r/r.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/octave/octave.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/ruby/ruby.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/sas/sas.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/javascript/javascript.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/shell/shell.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/julia/julia.js"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/codemirror/5.38.0/mode/markdown/markdown.js"></script>
      <script>
         {% include 'parts/sos-mode.js' %}
      </script>
      {{ doc.js() }}
      <script>
         CodeMirror.fromTextArea(document.getElementById("source-code"), {
           lineNumbers: true,
           styleActiveLine: true,
           matchBrackets: true,
           readOnly: true,
           theme: '{{ theme }}',
           mode: 'sos'
         });
         add_hoverdoc();
      </script>
   </body>
</html>
