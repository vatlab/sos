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
    <th> Worklow Name</th><td>{{workflow_name}}</td>
    </tr>
    <tr>
    <th> Worklow Start Time</th><td>{{workflow_start_time}}</td>
    </tr>
    <tr>
    <th> Worklow End Time</th><td>{{workflow_end_time}}</td>
    </tr>
    <tr>
    <th> Worklow Duration</th><td>{{workflow_duration}}</td>
    </tr>
    </table>

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

    {% for step, info in steps.items() %}
        <tr>
        <td>{{ info.stepname }}</td>
        <td>{{ info.input }}</td>
        <td>{{ info.output }}</td>
        <td>
            {% if info.completed.__substep_completed__ %} 
            {{ info.completed.__substep_completed__ }} completed &nbsp;
            {% endif %}
            {% if info.completed.__substep_skipped__ %} 
            {{ info.completed.__substep_skipped__ }} ignored
            {% endif %}
        </td>
        <td>
            {% if info.completed.__task_completed__ %} 
            {{ info.completed.__task_completed__ }} completed &nbsp;
            {% endif %}
            {% if info.completed.__task_skipped__ %} 
            {{ info.completed.__task_skipped__ }} ignored
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

    <h2> Summary </h2>
    <table class='summary_table'>
    <tr>
    {% if stat.__step_completed__  %}
    <tr>
    <th>Completed steps</th><td>{{ stat.__step_completed__ }}</td>
    </tr>
    {% endif %}
    {% if stat.__step_skipped__  %}
    <tr>
    <th>Ignored steps</th><td>{{ stat.__step_skipped__ }}</td>
    </tr>
    {% endif %}

    {% if stat.__substep_completed__  %}
    <tr>
    <th>Completed substeps</th><td>{{ stat.__substep_completed__ }}</td>
    </tr>
    {% endif %}
    {% if stat.__substep_skipped__  %}
    <tr>
    <th>Ignored substeps</th><td>{{ stat.__substep_skipped__ }}</td>
    </tr>
    {% endif %}

    {% if stat.__task_completed__  %}
    <tr>
    <th>Completed tasks</th><td>{{ stat.__task_completed__ }}</td>
    </tr>
    {% endif %}
    {% if stat.__task_skipped__  %}
    <tr>
    <th>Ignored tasks</th><td>{{ stat.__task_skipped__ }}</td>
    </tr>
    {% endif %}
    </table>
</body>
</html>


