var CodeCell = require('notebook/js/codecell').CodeCell;
var container_width=$('#site').width();
var my_scratchpad;


var Scratchpad = function (nb) {
    var scratchpad = this;
    this.notebook = nb;
    this.kernel = nb.kernel;
    this.km = nb.keyboard_manager;

    // create elements
    this.element = $("<div id='nbextension-scratchpad'>");
    // create my cell
    var cell = this.cell = new CodeCell(nb.kernel, {
      events: nb.events,
      config: nb.config,
      keyboard_manager: nb.keyboard_manager,
      notebook: nb,
      tooltip: nb.tooltip,
    });
    cell.set_input_prompt();
    this.element.append($("<div/>").addClass('cell-wrapper').append(this.cell.element));
    cell.render();
    cell.refresh();
    this.element.animate({
      height: 0,
    }, 100);
    this.cell.element.hide();


    // override ctrl/shift-enter to execute me if I'm focused instead of the notebook's cell
    var execute_and_select_action = this.km.actions.register({
      handler: $.proxy(this.execute_and_select_event, this),
    }, 'scratchpad-execute-and-select');
    var execute_action = this.km.actions.register({
      handler: $.proxy(this.execute_event, this),
    }, 'scratchpad-execute');
    var toggle_action = this.km.actions.register({
      handler: $.proxy(toggle_scratchpad, this),
    }, 'scratchpad-toggle');
    
    var shortcuts = {
      'shift-enter': execute_and_select_action,
      'ctrl-enter': execute_action,
      'ctrl-b': toggle_action,
    }
    this.km.edit_shortcuts.add_shortcuts(shortcuts);
    this.km.command_shortcuts.add_shortcuts(shortcuts);

    // finally, add me to the page
    $("body").append(this.element);
  };


  Scratchpad.prototype.execute_and_select_event = function (evt) {
    if (utils.is_focused(this.element)) {
      this.cell.execute();
    } else {
      this.notebook.execute_cell_and_select_below();
    }
  };

  Scratchpad.prototype.execute_event = function (evt) {
    if (utils.is_focused(this.element)) {
      this.cell.execute();
    } else {
      this.notebook.execute_selected_cells();
    }
  };


  function setup_scratchpad () {
    // lazy, hook it up to Jupyter.notebook as the handle on all the singletons
    console.log("Setting up scratchpad");
    return new Scratchpad(Jupyter.notebook);
  }


  function toggle_scratchpad(){
    if ($('#nbextension-scratchpad').height()==0){
        var site_height = $("#site").height();
        $('#nbextension-scratchpad').animate({
          height: site_height,
        }, 200);
        $('#nbextension-scratchpad .cell').show();
        $('#notebook-container').css('margin-left',$('#nbextension-scratchpad').width()+10);
        $('#notebook-container').css('margin-right',50);
        // $('#notebook-container').css('width',container_width-$('#nbextension-scratchpad').width()-30);
        $('#notebook-container').css('width','80%');
        $('.celltoolbar label').css('margin-left',0);
        $('.celltoolbar label').css('margin-right',0);
        my_scratchpad.cell.focus_editor()
    }else{
        $('#nbextension-scratchpad').animate({
              height: 0,
            }, 100);
        $('#nbextension-scratchpad .cell').hide();
        $('#notebook-container').css('margin-left','10px');
        $('#notebook-container').css('width',container_width);
        $('#notebook-container').css('margin-right',50);
    }

  }


  function load_scratchpad() {

    var load_css = function() {
        var css = document.createElement("style");
        css.type = "text/css";
        css.innerHTML = '#nbextension-scratchpad {position: absolute; left: 0; bottom: 0; width: 20%; background-color: #F8F5E1; border-left: 1px solid #aaa; border-top: 1px solid #aaa; z-index: 9000; } #nbextension-scratchpad .cell-wrapper {height: 100%; overflow: auto; } .scratchpad-btn {float: left; padding-right: 24px; opacity: 0.2; font-size: 24px; z-index: 9001; } .scratchpad-btn:hover {opacity: 1; } .scratchpad-close {display: none; position: absolute; float: right; bottom: 8px; right: 0; } .scratchpad-open {margin-top: -32px; }}'; 
        document.body.appendChild(css);
    };

    load_css();

    if (Jupyter.notebook.kernel) {
      my_scratchpad=setup_scratchpad();
    } else {
      events.on('kernel_ready.Kernel', setup_scratchpad);
    }
  }


function add_scratchpad_button() {
    if (!IPython.toolbar) {
      $([IPython.events]).on("app_initialized.NotebookApp", scratchpad_button);
      return;
    }
    if ($("#scratchpad_button").length === 0) {
      IPython.toolbar.add_buttons_group([
        {
          'label'   : 'Scratchpad',
          'icon'    : 'fa-cube',
          'callback':  toggle_scratchpad,
          'id'      : 'scratchpad_button'
        }
      ]);
    }
  };


function keepWidth(){
    if ($('#nbextension-scratchpad').height()!=0){
        $('#notebook-container').css('margin-left',$('#nbextension-scratchpad').width()+30);
        $('#notebook-container').css('width',container_width-$('#nbextension-scratchpad').width()-30);
        $('.celltoolbar label').css('margin-left',0);
        $('.celltoolbar label').css('margin-right',0);
    }
}

function patch_CodeCell_get_callbacks() {
    var previous_get_callbacks = CodeCell.prototype.get_callbacks;
    CodeCell.prototype.get_callbacks = function() {
        var that = this;
        var callbacks = previous_get_callbacks.apply(this, arguments);
        var prev_reply_callback = callbacks.shell.reply;
        callbacks.shell.reply = function(msg) {
            if (msg.msg_type === 'execute_reply') {
                keepWidth()
            }
            return prev_reply_callback(msg);
        };
        return callbacks;
    };
  }


