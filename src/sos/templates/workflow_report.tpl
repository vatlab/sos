<!DOCTYPE html>
<html lang="en">
   <head>
      <title>{{workflow_name}}</title>
      <style type="text/css">
         {% include "workflow_report.css" %}
      </style>
   </head>
   <body>
      <h1>SoS Workflow Report</h1>
      <table class='info_table'>
         {% if workflows[master_id].workflow_cmd %}
         <tr>
            <th> Command Line</th>
            <td><code>{{workflows[master_id].command_line}}</code></td>
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
      {% if subworkflows %}
      <h2>
      Subworkflows </h1>
      <table class='workflow_table'>
         <tr>
            <th>Workflow Name</th>
            <th>Start Time</th>
            <th>Duration</th>
            <th class="timeline-col">Timeline</th>
            <th>Steps</th>
            <th>Tasks</th>
         </tr>
         {% for name in subworkflows %}
         <tr>
            <td>{{ workflows[name].name }}</td>
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
               {% if workflows[name].stat.__task_completed__ %}
               {{ workflows[name].stat.__task_completed__  }} completed &nbsp;
               {% endif %}
               {% if workflows[name].stat.__task_skipped__  %}
               {{ workflows[name].stat.__task_skipped__ }} ignored
               {% endif %}
            </td>
         </tr>
         {% endfor %}
      </table>
      {% endif %}
      {% if steps %}
      <h2>
      Steps </h1>
      <table class='task_table'>
         <tr>
            <th>Step Name</th>
            <th>Input</th>
            <th>Output</th>
            <th>Start Time</th>
            <th>Duration</th>
            <th class="timeline-col">Timeline</th>
            <th>Substeps</th>
            <th>Tasks</th>
         </tr>
         {% for step in steps %}
         <tr>
            <td>{{ step.stepname }}</td>
            <td>{{ step.input }}</td>
            <td>{{ step.output }}</td>
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
      </table>
      {% endif %}
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
   </body>
</html>
