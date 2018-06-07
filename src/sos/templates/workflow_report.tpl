<!DOCTYPE html>
<html lang="en">
   <head>
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
      <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/css/bootstrap.min.css" integrity="sha384-WskhaSGFgHYWDcbwN70/dfYBj47jz9qbsMId/iRN3ewGhXQFZCSftd1LZCfmhktB" crossorigin="anonymous">

      <title>{{workflow_name}}</title>
      <style type="text/css">
         {% include "workflow_report.css" %}
      </style>
   </head>
   <body>
      <h1>SoS Workflow Report</h1>
      <table class='info_table'>
         {% if workflows[master_id].command_line %}
         <tr>
            <th> Command Line</th>
            <td><code>{{workflows[master_id].command_line}}</code></td>
         </tr>
         <tr>
            <th> Project directory</th>
            <td><code>{{workflows[master_id].project_dir}}</code></td>
         </tr>
         {% endif %}
         <tr>
            <th> Worklow Name</th>
            <td>{{ workflows[master_id].name }}</td>
         </tr>
         <tr>
            <th> Worklow Start Time</th>
            <td>{{ workflows[master_id].start_time_str }}</td>
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
      <h2> Steps </h2>
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
         <tr class="{{ "master_workflow" if name == master_id else "subworkflow" }}">
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
         <tr>
            <td></td>
            <td>{{ step.stepname }}</td>
            <td>{{ step.input_str }}</td>
            <td>{{ step.output_str }}</td>
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
      <h2> Tasks </h2>
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
            <td><code>{{ task }}</code></td>
            <td>{{ info.queue }}</td>
            <td>{{ "Ignored" if info.skipped else ("Success" if info.ret_code == 0 else "Failed") }}</td>
            <td><code>{{ info.tags }}</code></td>
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
      <h2> Execution DAG </h2>
      <img class="dag-image" src="data:image/png;base64,{{ workflows[master_id].dag_img }}">
      {% endif %}

    <footer>
     Report generated by <a href="https://vatlab.github.io/sos-docs/">SoS Workflow Engine</a> version {{sos_version}}, by {{ user }} on {{ time_now_str }}
    </footer>
    <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js" integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo" crossorigin="anonymous"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.3/umd/popper.min.js" integrity="sha384-ZMP7rVo3mIykV+2+9J3UJ46jBk0WLaUAdn689aCwoqbBJiSnjAK/l8WvCWPIPm49" crossorigin="anonymous"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.1/js/bootstrap.min.js" integrity="sha384-smHYKdLADwkXOn1EmN1qk/HfnUcbVRZyYmZ4qpPea6sjB/pTJ0euyQp0Mk8ck+5T" crossorigin="anonymous"></script>
    <script>
    $(document).ready(function(){
        $('[data-toggle="tooltip"]').tooltip();
    });
    </script>
 </body>
</html>
