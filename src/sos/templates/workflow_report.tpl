<!DOCTYPE html>
<html lang="en">
<head>
    <title>{{workflow_name}}</title>
</head>
<body>

    <h1>SoS Workflow Report</h1>

    <h2> Key information </h2>
    <ul>
    <li> Worklow Name: {{workflow_name}}</li>
    <li> Worklow Start Time: {{workflow_start_time}}</li>
    <li> Worklow End Time: {{workflow_end_time}}</li>
    <li> Worklow Duration: {{workflow_duration}}</li>
    </ul>

    {% if tasks %}
    <h2> Tasks </h1>
    <table class='task_table'>
    <tr>
    <th>Task ID</th>
    <th>Return Code</th>
    </tr>
    {% for task, info in tasks.items() %}
        <tr>
        <td>{{ task }}</td>
        <td>{{ info.ret_code }}</td>
        </tr>
    {% endfor %}
    </table>
    {% endif %}

    <h2> Summary </h2>
    <ul>
    {% if stat.__step_completed__  %}
    <li>Completed steps: {{ stat.__step_completed__ }}</li>
    {% endif %}
    {% if stat.__step_skipped__  %}
    <li>Ignored steps: {{ stat.__step_skipped__ }}</li>
    {% endif %}

    {% if stat.__substep_completed__  %}
    <li>Completed substeps: {{ stat.__substep_completed__ }}</li>
    {% endif %}
    {% if stat.__substep_skipped__  %}
    <li>Ignored substeps: {{ stat.__substep_skipped__ }}</li>
    {% endif %}

    {% if stat.__task_completed__  %}
    <li>Completed tasks: {{ stat.__task_completed__ }}</li>
    {% endif %}
    {% if stat.__task_skipped__  %}
    <li>Ignored tasks: {{ stat.__task_skipped__ }}</li>
    {% endif %}
    </ul>
</body>
</html>


