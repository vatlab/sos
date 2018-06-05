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
    {% if workflow_cmd %}
    <tr>
    <th> Command Line</th><td><code>{{workflow_cmd}}</code></td>
    </tr>
    {% endif %}
    <tr>
    <th> Worklow Name</th><td>{{ workflow_name[master_id] }}</td>
    </tr>
    <tr>
    <th> Worklow Start Time</th><td>{{ workflow_start_time[master_id] }}</td>
    </tr>
    <tr>
    <th> Worklow End Time</th><td>{{workflow_end_time[master_id]}}</td>
    </tr>
    <tr>
    <th> Worklow Duration</th><td>{{workflow_duration[master_id]}}</td>
    </tr>

    {% if workflow_stat[master_id].__step_completed__  %}
    <tr>
    <th>Completed steps</th><td>{{ workflow_stat[master_id].__step_completed__ }}</td>
    </tr>
    {% endif %}
    {% if workflow_stat[master_id].__step_skipped__  %}
    <tr>
    <th>Ignored steps</th><td>{{ workflow_stat[master_id].__step_skipped__ }}</td>
    </tr>
    {% endif %}

    {% if workflow_stat[master_id].__substep_completed__  %}
    <tr>
    <th>Completed substeps</th><td>{{ workflow_stat[master_id].__substep_completed__ }}</td>
    </tr>
    {% endif %}
    {% if workflow_stat[master_id].__substep_skipped__  %}
    <tr>
    <th>Ignored substeps</th><td>{{ workflow_stat[master_id].__substep_skipped__ }}</td>
    </tr>
    {% endif %}

    {% if workflow_stat[master_id].__task_completed__  %}
    <tr>
    <th>Completed tasks</th><td>{{ workflow_stat[master_id].__task_completed__ }}</td>
    </tr>
    {% endif %}
    {% if workflow_stat[master_id].__task_skipped__  %}
    <tr>
    <th>Ignored tasks</th><td>{{ workflow_stat[master_id].__task_skipped__ }}</td>
    </tr>
    {% endif %}

    </table>

    {% if subworkflows %}
    <h2> Subworkflows </h1>
    <table class='workflow_table'>
    <tr>
    <th>Workflow Name</th>
    <th>Start Time</th>
    <th>Duration</th>
    <th>Steps</th>
    <th>Tasks</th>
    </tr>

    {% for name in subworkflows %}
        <tr>
        <td>{{ workflow_name[name] }}</td>
        <td>{{ workflow_start_time[name] }}</td>
        <td>{{ workflow_duration[name] }}</td>
        <td>
            {% if workflow_stat[name].__step_completed__ %}
            {{ workflow_stat[name].__step_completed__  }} completed &nbsp;
            {% endif %}
            {% if workflow_stat[name].__step_skipped__  %}
            {{ workflow_stat[name].__step_skipped__ }} ignored
            {% endif %}
        </td>
        <td>
            {% if workflow_stat[name].__task_completed__ %}
            {{ workflow_stat[name].__task_completed__  }} completed &nbsp;
            {% endif %}
            {% if workflow_stat[name].__task_skipped__  %}
            {{ workflow_stat[name].__task_skipped__ }} ignored
            {% endif %}
        </td>
        </tr>
    {% endfor %}
    </table>
    {% endif %}

    {% if steps %}
    <h2> Steps </h1>
    <table class='task_table'>
    <tr>
    <th>Step Name</th>
    <th>Input</th>
    <th>Output</th>
    <th>Substeps</th>
    <th>Tasks</th>
    </tr>

    {% for step in steps %}
        <tr>
        <td>{{ step.stepname }}</td>
        <td>{{ step.input }}</td>
        <td>{{ step.output }}</td>
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
    <h2> Tasks </h1>
    <table class='task_table'>
    <tr>
    <th>Task ID</th>
    <th>Queue</th>
    <th>Status</th>
    <th>Tags</th>
    <th>Start Time</th>
    <th>Duration</th>
    <th>Peak CPU</th>
    <th>Peak RAM</th>
    </tr>
    {% for task, info in tasks.items() %}
        <tr>
        <td><code>{{ task }}</code></td>
        <td>{{ info.queue }}</td>
        <td>{{ "Ignored" if info.skipped else ("Success" if info.ret_code == 0 else "Failed") }}</td>
        <td><code>{{ info.tags }}</code></td>
        <td>{{ info.start_time }}</td>
        <td>{{ info.duration }}</td>
        <td>{{ info.peak_cpu }}</td>
        <td>{{ info.peak_mem }}</td>
        </tr>
    {% endfor %}
    </table>
    {% endif %}

    {% if dag_image %}
    <h2> Execution DAG </h2>

    <img class="dag-image" src="data:image/png;base64,{{ dag_image }}">

    {% endif %}

</body>
</html>
