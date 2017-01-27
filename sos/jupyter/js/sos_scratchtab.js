// This file is part of Script of Scripts (sos), a workflow system
// for the execution of commands and scripts in different languages.
// Please visit https://github.com/bpeng2000/SOS for more information.
//
// Copyright (C) 2016 Bo Peng (bpeng@mdanderson.org)
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
// This is a codemirror mode for SoS, currently a clone of Python codemirror.js
// with some minor modification. An expert on javascript and code mirror
// is required to make it work for other langauges that SoS supports.
//
// var CodeCell = require('notebook/js/codecell').CodeCell;
// var utils = require('base/js/utils');

var container_width = $('#site').width();
var my_scratchTab;

var scratchTab = function(nb) {
    var scratchTab = this;
    this.notebook = nb;
    this.kernel = nb.kernel;
    this.km = nb.keyboard_manager;

    // create elements
    this.element = $("<div id='nbextension-scratchTab'>");


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
    }, 'scratchTab-execute-and-select');
    var execute_action = this.km.actions.register({
        handler: $.proxy(this.execute_event, this),
    }, 'scratchTab-execute');
    var toggle_action = this.km.actions.register({
        handler: $.proxy(toggle_scratchTab, this),
    }, 'scratchTab-toggle');

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


scratchTab.prototype.execute_and_select_event = function(evt) {
    if (utils.is_focused(this.element)) {
        this.cell.execute();
    } else {
        this.notebook.execute_cell_and_select_below();
    }
};

scratchTab.prototype.execute_event = function(evt) {
    if (utils.is_focused(this.element)) {
        this.cell.execute();
    } else {
        this.notebook.execute_selected_cells();
    }
};


function setup_scratchTab() {
    // lazy, hook it up to Jupyter.notebook as the handle on all the singletons
    console.log("Setting up scratchTab");
    return new scratchTab(Jupyter.notebook);
}

function toggle_scratchTab() {
    if ($('#nbextension-scratchTab').height() < 1) {
        var site_height = $("#site").height();

        $('#nbextension-scratchTab').animate({
            height: site_height,
        }, 200);
        $('#nbextension-scratchTab .cell').show();
        $('#notebook-container').css('margin-left', $('#nbextension-scratchTab').width() + 10);
        $('#notebook-container').css('margin-right', 50);
        // $('#notebook-container').css('width',container_width-$('#nbextension-scratchTab').width()-30);
        $('#notebook-container').css('width', '80%');
        $('.celltoolbar label').css('margin-left', 0);
        $('.celltoolbar label').css('margin-right', 0);
        my_scratchTab.cell.focus_editor()
    } else {
        $('#nbextension-scratchTab').animate({
            height: 0,
        }, 100);
        $('#nbextension-scratchTab .cell').hide();
        $('#notebook-container').css('margin-left', '10px');
        $('#notebook-container').css('width', container_width - 20);
        $('#notebook-container').css('margin-right', 50);
    }
}

function load_scratchTab() {

    var load_css = function() {
        var css = document.createElement("style");
        css.type = "text/css";
        css.innerHTML = '#nbextension-scratchTab {position: absolute; left: 0; bottom: 0; width: 20%; background-color: #F8F5E1; border-left: 1px solid #aaa; border-top: 1px solid #aaa; z-index: 9000; } #nbextension-scratchTab .cell-wrapper {height: 100%; overflow: auto; } .scratchTab-btn {float: left; padding-right: 24px; opacity: 0.2; font-size: 24px; z-index: 9001; } .scratchTab-btn:hover {opacity: 1; } .scratchTab-close {display: none; position: absolute; float: right; bottom: 8px; right: 0; } .scratchTab-open {margin-top: -32px; }}';
        document.body.appendChild(css);
    };

    load_css();

    if (Jupyter.notebook.kernel) {
        my_scratchTab = setup_scratchTab();
    } else {
        events.on('kernel_ready.Kernel', setup_scratchTab);
    }
}

function add_scratchTab_button() {
    if (!IPython.toolbar) {
        $([IPython.events]).on("app_initialized.NotebookApp", scratchTab_button);
        return;
    }
    if ($("#scratchTab_button").length === 0) {
        IPython.toolbar.add_buttons_group([{
            'label': 'scratchTab',
            'icon': 'fa-cube',
            'callback': toggle_scratchTab,
            'id': 'scratchTab_button'
        }]);
    }
};

function keepWidth() {
    if ($('#nbextension-scratchTab').height() != 0) {
        $('#notebook-container').css('margin-left', $('#nbextension-scratchTab').width() + 30);
        $('#notebook-container').css('width', container_width - $('#nbextension-scratchTab').width() - 30);
        $('.celltoolbar label').css('margin-left', 0);
        $('.celltoolbar label').css('margin-right', 0);
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
