{% import 'parts/hover_doc.tpl' as doc %}

<!DOCTYPE html>
<html lang="en">
   <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
      <link type="text/css" rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/jquery-dropdown/2.0.3/jquery.dropdown.css" />
      <link rel="stylesheet" href="https://codemirror.net/lib/codemirror.css">
      <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/css/bootstrap.min.css" integrity="sha384-WskhaSGFgHYWDcbwN70/dfYBj47jz9qbsMId/iRN3ewGhXQFZCSftd1LZCfmhktB" crossorigin="anonymous">
      <title>{{workflow_name}}</title>
      <style type="text/css">
         {% include "parts/sos_report.css" %}
      </style>
   </head>
   <body>
      <h1 class='mt-0'>SoS Workflow Report</h1>
      <table class='info_table'>
         {% if workflows[master_id].command_line %}
         <tr>
            <th> Command Line</th>
            <td><kbd>{{workflows[master_id].command_line}}</kbd></td>
         </tr>
         <tr>
            <th> Project directory</th>
            <td><samp>{{workflows[master_id].project_dir}}</samp></td>
         </tr>
         {% endif %}
         <tr>
            <th> Worklow Name</th>
            <td><samp>{{ workflows[master_id].name }}</samp></td>
         </tr>
         <tr>
            <th> Worklow ID</th>
            <td><samp>{{ master_id }}</samp></td>
         </tr>
         <tr>
            <th> Worklow Start Time</th>
            <td>{{ workflows[master_id].start_time_str }}
               (<span id="workflow-start-time">{{ workflows[master_id].start_time }}</span>)
            </td>
         </tr>
         <tr>
            <th> Worklow End Time</th>
            <td>{{ workflows[master_id].end_time_str }}</td>
         </tr>
         <tr>
            <th> Worklow Duration</th>
            <td>{{ workflows[master_id].duration_str }}</td>
         </tr>
         {% if workflows[master_id].stat.__step_completed__  %}
         <tr>
            <th>Completed steps</th>
            <td>{{ workflows[master_id].stat.__step_completed__ }}</td>
         </tr>
         {% endif %}
         {% if workflows[master_id].stat.__step_skipped__  %}
         <tr>
            <th>Ignored steps</th>
            <td>{{ workflows[master_id].stat.__step_skipped__ }}</td>
         </tr>
         {% endif %}
         {% if workflows[master_id].stat.__substep_completed__  %}
         <tr>
            <th>Completed substeps</th>
            <td>{{ workflows[master_id].stat.__substep_completed__ }}</td>
         </tr>
         {% endif %}
         {% if workflows[master_id].stat.__substep_skipped__  %}
         <tr>
            <th>Ignored substeps</th>
            <td>{{ workflows[master_id].stat.__substep_skipped__ }}</td>
         </tr>
         {% endif %}
         {% if workflows[master_id].stat.__task_completed__  %}
         <tr>
            <th>Completed tasks</th>
            <td>{{ workflows[master_id].stat.__task_completed__ }}</td>
         </tr>
         {% endif %}
         {% if workflows[master_id].stat.__task_skipped__  %}
         <tr>
            <th>Ignored tasks</th>
            <td>{{ workflows[master_id].stat.__task_skipped__ }}</td>
         </tr>
         {% endif %}
      </table>
      <h2 class='mt-3'> Steps </h2>
      <table class='Steps'>
         <tr>
            <th>Workflow</th>
            <th>Step</th>
            <th>Input</th>
            <th>Output</th>
            <th>Start Time</th>
            <th>Duration</th>
            <th class="timeline-col">Timeline</th>
            <th>Steps</th>
            <th>Substeps</th>
            <th>Tasks</th>
         </tr>
         {% for name in workflows.keys() %}
         <tr class="{{ "master_workflow" if name == master_id else "subworkflow"}} workflow-row">
         <th>{{ workflows[name].name }}</th>
         <td></td>
         <td></td>
         <td></td>
         <td>{{ workflows[name].start_time_str }}</td>
         <td>{{ workflows[name].duration_str }}</td>
         <td>
            <div class="timeline-cell">
               <div class="timeline-before" style="width:{{ workflows[name].before_percent }}%"></div>
               <div class="timeline-during" style="width:{{ workflows[name].during_percent }}%"></div>
               <div class="timeline-after" style="width:{{ workflows[name].after_percent }}%"></div>
            </div>
         </td>
         <td>
            {% if workflows[name].stat.__step_completed__ %}
            {{ workflows[name].stat.__step_completed__  }} completed &nbsp;
            {% endif %}
            {% if workflows[name].stat.__step_skipped__  %}
            {{ workflows[name].stat.__step_skipped__ }} ignored
            {% endif %}
         </td>
         <td>
            {% if workflows[name].stat.__substep_completed__ %}
            {{ workflows[name].stat.__substep_completed__  }} completed &nbsp;
            {% endif %}
            {% if workflows[name].stat.__substep_skipped__  %}
            {{ workflows[name].stat.__substep_skipped__ }} ignored
            {% endif %}
         </td>
         <td>
            {% if workflows[name].stat.__task_completed__ %}
            {{ workflows[name].stat.__task_completed__  }} completed &nbsp;
            {% endif %}
            {% if workflows[name].stat.__task_skipped__  %}
            {{ workflows[name].stat.__task_skipped__ }} ignored
            {% endif %}
         </td>
         </tr>
         {% for step in steps[name] %}
         <tr class="{{ "master_workflow" if name == master_id else "subworkflow"}} step-row">
         <td></td>
         <td>{{ step.stepname }}</td>
         <td>
            {% if step.input %}
            <div class="dropdown">
               <span data-toggle="dropdown" class="filelist">
               {{ step.input[0][0] | basename | truncate(20) }}{% if step.input|length > 1 %} ({{ step.input|length }}) {% endif %}
               </span>
               <ul class="dropdown-menu">
                  {% for file in step.input %}
                  <li> {{file[0]}} <span class="text-muted filesize">{{file[1]|filesizeformat}}</span>
                  <li>
                     {% endfor %}
               </ul>
            </div>
            {% endif %}
         </td>
         <td>
            {% if step.output %}
            <div class="dropdown">
               <span data-toggle="dropdown" class="filelist">
               {{ step.output[0][0] | basename | truncate(20) }}{% if step.output|length > 1 %} ({{ step.output|length }}) {% endif %}
               </span>
               <ul class="dropdown-menu">
                  {% for file in step.output %}
                  <li> {{file[0]}} <span class="text-muted filesize">{{file[1] | filesizeformat}}</span>
                  <li>
                     {% endfor %}
               </ul>
            </div>
            {% endif %}
         </td>
         <td>{{ step.start_time_str }}</td>
         <td>{{ step.duration_str }}</td>
         <td>
            <div class="timeline-cell">
               <div class="timeline-before" style="width:{{ step.before_percent }}%"></div>
               <div class="timeline-during" style="width:{{ step.during_percent }}%"></div>
               <div class="timeline-after" style="width:{{ step.after_percent }}%"></div>
            </div>
         </td>
         <td>
            {% if step.completed.__step_completed__ %}
            {{ step.completed.__step_completed__ }} completed &nbsp;
            {% endif %}
            {% if step.completed.__step_skipped__ %}
            {{ step.completed.__step_skipped__ }} ignored
            {% endif %}
         </td>
         <td>
            {% if step.completed.__substep_completed__ %}
            {{ step.completed.__substep_completed__ }} completed &nbsp;
            {% endif %}
            {% if step.completed.__substep_skipped__ %}
            {{ step.completed.__substep_skipped__ }} ignored
            {% endif %}
         </td>
         <td>
            {% if step.completed.__task_completed__ %}
            {{ step.completed.__task_completed__ }} completed &nbsp;
            {% endif %}
            {% if step.completed.__task_skipped__ %}
            {{ step.completed.__task_skipped__ }} ignored
            {% endif %}
         </td>
         </tr>
         {% endfor %}
         {% endfor %}
      </table>
      {% if tasks %}
      <h2 class='mt-3'> Tasks </h2>
      <table class='task_table'>
         <tr>
            <th>Task ID</th>
            <th>Queue</th>
            <th>Status</th>
            <th>Tags</th>
            <th>Start Time</th>
            <th>Duration</th>
            <th class="timeline-col">Timeline</th>
            <th>Peak CPU</th>
            <th>Peak RAM</th>
         </tr>
         {% for task, info in tasks.items() %}
         <tr>
            <td><samp>{{ task }}</samp></td>
            <td>{{ info.queue }}</td>
            <td>{{ "Ignored" if info.skipped else ("Success" if info.ret_code == 0 else "Failed") }}</td>
            <td><samp>{{ info.tags }}</samp></td>
            <td>{{ info.start_time_str }}</td>
            <td>{{ info.duration_str }}</td>
            <td>
               <div class="timeline-cell">
                  <div class="timeline-before" style="width:{{ info.before_percent }}%"></div>
                  <div class="timeline-during" style="width:{{ info.during_percent }}%"></div>
                  <div class="timeline-after" style="width:{{ info.after_percent }}%"></div>
               </div>
            </td>
            <td>{{ info.peak_cpu_str }}</td>
            <td>{{ info.peak_mem_str }}</td>
         </tr>
         {% endfor %}
      </table>
      {% endif %}
      {% if workflows[master_id].dag_img %}
      <h2 class='mt-3'> Execution DAG </h2>
      <img class="dag-image" src="data:image/png;base64,{{ workflows[master_id].dag_img }}">
      {% endif %}
      <h2 class='mt-3'> Script </h2>
      <div class="file">
         <div class="fileheader">
            <div class="fileinfo">
               {{ workflows[master_id].script.splitlines() | length }} lines
               <span class="file-info-divider"></span>
               {{ workflows[master_id].script | length | filesizeformat}}
            </div>
         </div>
         <div class="filecontent">
            <textarea rows="{{ workflows[master_id].script.splitlines() | count }}"
			id="source-code" class="sos-source" name="code">{{ workflows[master_id].script }}</textarea>
         </div>
      </div>
      {% if transcripts %}
      <h2 class='mt-3'> Transcript </h2>
      <div class="transcript-list">
        {% for id, trans in transcripts.items() %}
        <div class="transcript-stepname">{{ id }}</div>
        {% for tran in trans %}
        <div class="transcript">
            <span class="transcript-time">{{ tran.start_time_str }}</span>
            <div class="transcript-command">{{ tran.command }}</div>
            <div class="transcript-script"><pre>{{ tran.script | trim }}</pre></div>
        </div>
        {% endfor %}
        {% endfor %}
      </div>
      {% endif %}
      <footer>
         <a class="sos-logo" href="https://vatlab.github.io/sos-docs">
         <img src="http://vatlab.github.io/sos-docs/img/sos_icon.svg" alt="sos_icon">
         </a>
         &nbsp;
         Report generated by <a href="https://vatlab.github.io/sos-docs/">SoS Workflow Engine</a> version <samp>{{sos_version}}</samp>, by <samp>{{ user }}</samp> on {{ time_now_str }}
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
         function timeElapsed(start_date)
         {
           let now = new Date()
           let seconds = parseInt(now / 1000 - start_date);

           let year = now.getUTCFullYear() - new Date(start_date * 1000).getUTCFullYear();
           if (year > 0) {
             return year + " year" + (year > 1 ? "s" : "") + " ago";
           }
           var day = Math.floor(seconds / 86400);
           if (day > 0) {
             return day + " day" + (day > 1 ? "s" : "") + " ago";
           }
           var hh = Math.floor((seconds % 86400) / 3600);
           if (hh > 0) {
             return hh + " hour" + (hh > 1 ? "s" : "") + " ago";
           }
           var mm = Math.floor((seconds % 3600) / 60);
           if (mm > 0) {
             return mm + " minute" + (mm > 1 ? "s" : "") + " ago";
           }
           return "less than 1 minute ago"
         };
         $(document).ready(function(){
           $('[data-toggle="tooltip"]').tooltip();
           $(".workflow-row").click(function(){
             $(this).nextUntil('tr.workflow-row').slideToggle(100, function(){
            });
          });
          $("#workflow-start-time").text(timeElapsed($("#workflow-start-time")[0].innerText));
         });

         CodeMirror.fromTextArea(document.getElementById("source-code"), {
           lineNumbers: true,
           styleActiveLine: true,
           matchBrackets: true,
           readOnly: true,
           indentUnit: 5,
           mode: 'sos'
         });
         add_hoverdoc();
      </script>
   </body>
</html>
