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
define(['jquery'], function($) {

    "use strict";
    //variables defined as global which enable access from imported scripts.
    window.BackgroundColor = {}
    window.DisplayName = {}
    window.KernelName = {}
    window.KernelList = []
    window.events = require('base/js/events');
    window.utils = require('base/js/utils');
    window.Jupyter = require('base/js/namespace');
    window.CodeCell = require('notebook/js/codecell').CodeCell;

    window.default_kernel = 'sos';
    window.kernel_updated = false;

    // initialize BackgroundColor etc from cell meta data
    if (!('sos' in IPython.notebook.metadata))
        IPython.notebook.metadata['sos'] = {
            'kernels': [
                ['sos', 'SoS', '']
            ]
        };

    window.data = IPython.notebook.metadata['sos']['kernels'];
    for (var i = 0; i < data.length; i++) {
        // BackgroundColor is color
        BackgroundColor[data[i][0]] = data[i][2];
        BackgroundColor[data[i][1]] = data[i][2];
        // DisplayName
        DisplayName[data[i][0]] = data[i][1];
        DisplayName[data[i][1]] = data[i][1];
        // Name
        KernelName[data[i][0]] = data[i][0];
        KernelName[data[i][1]] = data[i][0];
        // KernelList, use displayed name
        KernelList.push([data[i][1], data[i][1]])
    }

    var my_execute = function(code, callbacks, options) {
        "use strict"
        var cells = IPython.notebook.get_cells();
        for (var i = cells.length - 1; i >= 0; --i) {
            // this is the cell that is being executed...
            // according to this.set_input_prompt('*') before execute is called.
            // also, because a cell might be starting without a previous cell
            // being finished, we should start from reverse and check actual code
            if (cells[i].input_prompt_number == '*' && code == cells[i].get_text()) {
                // use cell kernel if meta exists, otherwise use window.default_kernel
                return this.orig_execute(
                    // passing to kernel
                    // 1. the default kernel (might have been changed from menu bar
                    // 2. cell kernel (might be unspecified for new cell)
                    // 3. cell index (for setting style after execution)
                    // in addition, the frontend command will send a "--list-kernel" request if
                    // the frontend is not correctly initialized, possibly because the kernel was
                    // not ready when the frontend sent the command `%listkernel`.
                    "%frontend " +
                    (window.kernel_updated ? "" : " --list-kernel ") +
                    (window.my_panel.displayed ? " --use-panel" : "") +
                    " --default-kernel " + window.default_kernel +
                    " --cell-kernel " + cells[i].metadata.kernel +
                    " --cell " + i.toString() + "\n" + code,
                    callbacks, options)
            }
        }
        // if this is a command from scratch pad (not part of the notebook)
        return this.orig_execute(
            "%frontend " +
            (window.kernel_updated ? "" : " --list-kernel ") +
            " --default-kernel " + window.default_kernel +
            " --cell-kernel " + window.my_panel.cell.metadata.kernel +
            " --cell -1 " + "\n" + code,
            callbacks, {
                'silent': false,
                'store_history': false
            })
    }

    function register_sos_comm() {
        // comm message sent from the kernel
        Jupyter.notebook.kernel.comm_manager.register_target('sos_comm',
            function(comm, msg) {
                comm.on_msg(function(msg) {
                    // when the notebook starts it should receive a message in the format of
                    // a nested array of elements such as
                    //
                    // "ir", "R", "#ABackgroundColorDEF"
                    //
                    // where are kernel name (jupyter kernel), displayed name (SoS), and background
                    // color assigned by the language module. The user might use name ir or R (both
                    // acceptable) but the frontend should only display displayed name, and send
                    // the real kernel name back to kernel (%frontend and metadata).
                    //
                    // there are two kinds of messages from my_execute
                    // 1. cell_idx: kernel
                    //     the kernel used for the cell with source
                    // 2. None: kernel
                    //     the kernel for the new cell

                    var data = msg.content.data;
                    var msg_type = msg.metadata.msg_type;

                    if (msg_type == 'kernel-list') {
                        if (window.kernel_updated)
                            return;
                        for (var i = 0; i < data.length; i++) {
                            // BackgroundColor is color
                            BackgroundColor[data[i][0]] = data[i][2];
                            BackgroundColor[data[i][1]] = data[i][2];
                            // DisplayName
                            DisplayName[data[i][0]] = data[i][1];
                            DisplayName[data[i][1]] = data[i][1];
                            // Name
                            KernelName[data[i][0]] = data[i][0];
                            KernelName[data[i][1]] = data[i][0];
                            // KernelList, use displayed name
                            if (KernelList.findIndex((item) => item[0] === data[i][1]) == -1)
                                KernelList.push([data[i][1], data[i][1]]);
                            // if the kernel is not in metadata, push it in
                            var k_idx = IPython.notebook.metadata['sos']['kernels'].findIndex((item) => item[0] === data[i][0])
                            if (k_idx == -1)
                                IPython.notebook.metadata['sos']['kernels'].push(data[i])
                            else {
                                // if language exist update the display name and color, in case it was using old ones
                                IPython.notebook.metadata['sos']['kernels'][k_idx][1] = data[i][1];
                                IPython.notebook.metadata['sos']['kernels'][k_idx][2] = data[i][2];
                            }
                        }
                        //add dropdown menu of kernels in frontend
                        load_select_kernel();
                        window.kernel_updated = true;
                        console.log('kernel list updated');
                    } else if (msg_type == 'default-kernel') {
                        // update the cells when the notebook is being opened.
                        // we also set a global kernel to be used for new cells
                        $('#kernel_selector').val(DisplayName[data]);
                        // a side effect of change is cells without metadata kernel info will change background
                        $('#kernel_selector').change();
                    } else if (msg_type == 'cell-kernel') {
                        // get cell from passed cell index, which was sent through the
                        // %frontend magic
                        if (data[0] == -1)
                            var cell = window.my_panel.cell;
                        else
                            var cell = IPython.notebook.get_cell(data[0]);
                        if (cell.metadata.kernel != KernelName[data[1]]) {
                            cell.metadata.kernel = KernelName[data[1]];
                            // set meta information
                            changeStyleOnKernel(cell, data[1])
                        }
                    } else if (msg_type == 'preview-input') {
                        cell = window.my_panel.cell;
                        cell.clear_input();
                        cell.set_text(data);
                        cell.clear_output();
                    } else if (msg_type == 'preview-kernel') {
                        changeStyleOnKernel(window.my_panel.cell, data);
                    } else {
                        // this is preview output
                        cell = window.my_panel.cell;
                        data.output_type = msg_type;
                        cell.output_area.append_output(data);
                        /*
                        if (msg_type === 'display_data')
                            cell.output_area.append_display_data(data);
                        else if (msg_type === 'stream')
                            cell.output_area.append_stream(data);
                        else
                            cell.output_area.append_stream('Unknown msg type ' + msg_type + '\nPlease notify maintainer of SoS of this bug.');
                        */
                        // remove output prompt
                        while (true) {
                            var op = cell.element[0].getElementsByClassName('out_prompt_overlay');
                            if (op.length > 0)
                                op[0].parentNode.removeChild(op[0]);
                            else
                                break;
                        }
                        while (true) {
                            var op = cell.element[0].getElementsByClassName('prompt');
                            if (op.length > 0)
                                op[0].parentNode.removeChild(op[0]);
                            else
                                break;
                        }
                        var ops = cell.element[0].getElementsByClassName('output_subarea');
                        for (var op = 0; op < ops.length; op++)
                            ops[op].style.maxWidth = '100%';
                        cell.output_area.expand();
                    }
                });
            }
        );
        console.log('sos comm registered');
    }


    function wrap_execute() {
        if (!window.kernel_updated) {
            IPython.notebook.kernel.execute('%frontend --list-kernel', [], {
                'silent': true,
                'store_history': false
            });
            console.log('kernel list requested');
        }
        // override kernel execute with the wrapper.
        IPython.notebook.kernel.orig_execute = IPython.notebook.kernel.execute;
        IPython.notebook.kernel.execute = my_execute;
    }

    function changeStyleOnKernel(cell, type) {
        var sel = cell.element[0].getElementsByTagName('select')[0]
        var opts = sel.options;
        for (var opt, j = 0; opt = opts[j]; j++) {
            if (opt.value == DisplayName[type]) {
                sel.selectedIndex = j;
                break;
            }
        }
        var ip = cell.element[0].getElementsByClassName('input_prompt');
        var op = cell.element[0].getElementsByClassName('out_prompt_overlay');

        // cell in panel does not have prompt area
        if (ip.length == 0) {
            if (BackgroundColor[type])
                cell.element[0].getElementsByClassName('input')[0].style.backgroundColor = BackgroundColor[type];
            else
                cell.element[0].getElementsByClassName('input')[0].style.backgroundColor = '';
            return;
        }

        if (BackgroundColor[type]) {
            ip[0].style.backgroundColor = BackgroundColor[type];
            op[0].style.backgroundColor = BackgroundColor[type];
        } else {
            // Use '' to remove background-color?
            ip[0].style.backgroundColor = '';
            op[0].style.backgroundColor = '';
        }

    }


    function load_select_kernel() {
        // this function will be called twice, the first time when the notebook is loaded
        // to create UT elements using the information from notebook metadata. The second
        // time will be caused whent the backend sends the frontend a list of available kernels
        // this is why we should not add additional UI elements when the function is called
        // the second time.

        //change css for CellToolBar
        var load_css = function() {
            var css = document.createElement("style");
            css.type = "text/css";
            css.innerHTML = '.code_cell .celltoolbar {width:10%;background:none;border:none;border-bottom:none;z-index: 1000;position:relative;margin-bottom:-50pt;float:right;}  .text_cell .celltoolbar {display:none}';
            document.body.appendChild(css);
        };

        load_css();

        var CellToolbar = IPython.CellToolbar;
        // the cell tool bar might have been added by the previous load_select_kernel call
        var slideshow_preset = [];
        var select_type = CellToolbar.utils.select_ui_generator(
            KernelList,
            // setter
            function(cell, value) {
                // we check that the slideshow namespace exist and create it if needed
                //if (cell.metadata.kernel == undefined) {
                cell.metadata.kernel = KernelName[value];
                var ip = cell.element[0].getElementsByClassName('input_prompt');
                var op = cell.element[0].getElementsByClassName('out_prompt_overlay');
                // cell in panel does not have prompt area
                if (ip.length == 0) {
                    if (BackgroundColor[value])
                        cell.element[0].getElementsByClassName('input')[0].style.backgroundColor = BackgroundColor[value];
                    else
                        cell.element[0].getElementsByClassName('input')[0].style.backgroundColor = '';
                    return;
                }

                //}
                // cell.metadata.kernel = KernelName[value];
                if (BackgroundColor[value]) {
                    ip[0].style.backgroundColor = BackgroundColor[value];
                    op[0].style.backgroundColor = BackgroundColor[value];
                } else {
                    // Use '' to remove background-color?
                    ip[0].style.backgroundColor = '';
                    op[0].style.backgroundColor = '';
                }
            },
            //geter
            function(cell) {
                var ns = cell.metadata.kernel;
                return (ns == undefined) ? undefined : ns.kernel
            },
            "");

        if (CellToolbar.list_presets().indexOf("Select cell kernel") > 0)
            CellToolbar.unregister_preset('Select cell kernel');
        CellToolbar.register_callback('slideshow.select', select_type);
        slideshow_preset.push('slideshow.select');
        var reveal_preset = slideshow_preset.slice();
        CellToolbar.register_preset('Select cell kernel', reveal_preset);
        // console.log('Select cell kernel loaded.');
        CellToolbar.global_show();
        CellToolbar.activate_preset('Select cell kernel');

        var cells = IPython.notebook.get_cells();
        for (var i in cells) {
            if (cells[i].cell_type == 'code') {
                changeStyleOnKernel(cells[i], cells[i].metadata.kernel);
            }
        }

        var dropdown = $("<select></select>").attr("id", "kernel_selector")
            .css("margin-left", "0.75em")
            .attr("class", "form-control select-xs")
        // .change(select_kernel);
        if (Jupyter.toolbar.element.has('#kernel_selector').length == 0)
            Jupyter.toolbar.element.append(dropdown);
        // remove any existing items
        $('#kernel_selector').empty();
        $.each(KernelList, function(key, value) {
            $('#kernel_selector')
                .append($("<option></option>")
                    .attr("value", DisplayName[value[0]])
                    .text(DisplayName[value[0]]));
        });
        $('#kernel_selector').val("SoS");
        $('#kernel_selector').change(function() {
            var kernel_type = $("#kernel_selector").val();

            window.default_kernel = kernel_type;

            var cells = IPython.notebook.get_cells();
            for (var i in cells) {
                if (cells[i].cell_type == 'code' && !cells[i].metadata.kernel) {
                    changeStyleOnKernel(cells[i], kernel_type);
                }
            }
        });
    }

    function changeCellStyle() {
        var cells = IPython.notebook.get_cells();
        // setting up background color and selection according to notebook metadata
        for (var i in cells) {
            if (cells[i].cell_type == 'code') {
                changeStyleOnKernel(cells[i], cells[i].metadata.kernel);
            }
        }
    }

    var container_width = $('#site').width();
    var my_panel;

    var panel = function(nb) {
        var panel = this;
        this.notebook = nb;
        this.kernel = nb.kernel;
        this.km = nb.keyboard_manager;

        // create elements
        this.element = $("<div id='sos-panel'>");

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

        // remove input and output prompt to save some space
        var ip = this.cell.element[0].getElementsByClassName('input_prompt');
        ip[0].parentNode.removeChild(ip[0]);
        // move the language selection stuff to the top
        this.cell.element[0].getElementsByClassName('celltoolbar')[0].style.marginBottom = 0;

        // override ctrl/shift-enter to execute me if I'm focused instead of the notebook's cell
        var execute_and_select_action = this.km.actions.register({
            handler: $.proxy(this.execute_and_select_event, this),
        }, 'panel-execute-and-select');
        var execute_action = this.km.actions.register({
            handler: $.proxy(this.execute_event, this),
        }, 'panel-execute');
        var toggle_action = this.km.actions.register({
            handler: $.proxy(toggle_panel, this),
        }, 'panel-toggle');

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


    panel.prototype.execute_and_select_event = function(evt) {
        if (utils.is_focused(this.element)) {
            this.cell.execute();
        } else {
            this.notebook.execute_cell_and_select_below();
        }
    };

    panel.prototype.execute_event = function(evt) {
        if (utils.is_focused(this.element)) {
            this.cell.execute();
        } else {
            this.notebook.execute_selected_cells();
        }
    };


    function setup_panel() {
        // lazy, hook it up to Jupyter.notebook as the handle on all the singletons
        console.log("Setting up panel");
        return new panel(Jupyter.notebook);
    }

    function toggle_panel() {
        if ($('#sos-panel').height() < 1) {
            var site_height = $("#site").height();

            $('#sos-panel').animate({
                height: site_height,
            }, 200);
            $('#sos-panel .cell').show();
            $('#notebook-container').css('margin-left', $('#sos-panel').width() + 10);
            $('#notebook-container').css('margin-right', 50);
            // $('#notebook-container').css('width',container_width-$('#sos-panel').width()-30);
            $('#notebook-container').css('width', '80%');
            $('.celltoolbar label').css('margin-left', 0);
            $('.celltoolbar label').css('margin-right', 0);
            my_panel.cell.focus_editor()
            my_panel.displayed = true;
        } else {
            $('#sos-panel').animate({
                height: 0,
            }, 100);
            $('#sos-panel .cell').hide();
            $('#notebook-container').css('margin-left', '10px');
            $('#notebook-container').css('width', container_width - 20);
            $('#notebook-container').css('margin-right', 50);
            my_panel.displayed = false;
        }
    }

    function load_panel() {

        var load_css = function() {
            var css = document.createElement("style");
            css.type = "text/css";
            css.innerHTML = '#sos-panel {position: absolute; left: 0; bottom: 0; width: 25%; background-color: #F8F5E1; border-left: 1px solid #aaa; border-top: 1px solid #aaa; z-index: 9000; } #sos-panel .cell-wrapper {height: 100%; overflow: auto; } .panel-btn {float: left; padding-right: 24px; opacity: 0.2; font-size: 24px; z-index: 9001; } .panel-btn:hover {opacity: 1; } .panel-close {display: none; position: absolute; float: right; bottom: 8px; right: 0; } .panel-open {margin-top: -32px; }}';
            document.body.appendChild(css);
        };

        load_css();

        if (Jupyter.notebook.kernel) {
            my_panel = setup_panel();
        } else {
            events.on('kernel_ready.Kernel', setup_panel);
        }
    }

    function add_panel_button() {
        if (!IPython.toolbar) {
            $([IPython.events]).on("app_initialized.NotebookApp", panel_button);
            return;
        }
        if ($("#panel_button").length === 0) {
            IPython.toolbar.add_buttons_group([{
                'label': 'scratch tab',
                'icon': 'fa-cube',
                'callback': toggle_panel,
                'id': 'panel_button'
            }]);
        }
    };

    function adjustPanel() {
        if ($('#sos-panel').height() != 0) {
            $('#notebook-container').css('margin-left', $('#sos-panel').width() + 30);
            $('#notebook-container').css('width', container_width - $('#sos-panel').width() - 30);
            $('.celltoolbar label').css('margin-left', 0);
            $('.celltoolbar label').css('margin-right', 0);
        }
        // remove output prompt of the cell panel
        cell = window.my_panel.cell;
        // remove output prompt
        while (true) {
            var op = cell.element[0].getElementsByClassName('out_prompt_overlay');
            if (op.length > 0)
                op[0].parentNode.removeChild(op[0]);
            else
                break;
        }
        while (true) {
            var op = cell.element[0].getElementsByClassName('prompt');
            if (op.length > 0)
                op[0].parentNode.removeChild(op[0]);
            else
                break;
        }
        var ops = cell.element[0].getElementsByClassName('output_subarea');
        for (var op = 0; op < ops.length; op++)
            ops[op].style.maxWidth = '100%';
        cell.output_area.expand();
    }

    function patch_CodeCell_get_callbacks() {
        var previous_get_callbacks = CodeCell.prototype.get_callbacks;
        CodeCell.prototype.get_callbacks = function() {
            var that = this;
            var callbacks = previous_get_callbacks.apply(this, arguments);
            var prev_reply_callback = callbacks.shell.reply;
            callbacks.shell.reply = function(msg) {
                if (msg.msg_type === 'execute_reply') {
                    adjustPanel()
                }
                return prev_reply_callback(msg);
            };
            return callbacks;
        };
    }

    var onload = function() {

        // setting up frontend using existing metadata (without executing anything)
        load_select_kernel();
        changeCellStyle();
        events.on('kernel_connected.Kernel', register_sos_comm);
        events.on('kernel_ready.Kernel', wrap_execute);

        load_panel();
        add_panel_button();
        patch_CodeCell_get_callbacks();

        // define SOS CodeMirror syntax highlighter
        (function(mod) {
            //if (typeof exports == "object" && typeof module == "object") // CommonJS
            // mod(require("../../lib/codemirror"));
            //else if (typeof define == "function" && define.amd) // AMD
            //  define(["../../lib/codemirror"], mod);
            //else // Plain browser env
            mod(CodeMirror);
        })(function(CodeMirror) {
            "use strict";

            function wordRegexp(words) {
                return new RegExp("^((" + words.join(")|(") + "))\\b");
            }

            var wordOperators = wordRegexp(["and", "or", "not", "is"]);
            var commonKeywords = ["as", "assert", "break", "class", "continue",
                "def", "del", "elif", "else", "except", "finally",
                "for", "from", "global", "if", "import",
                "lambda", "pass", "raise", "return",
                "try", "while", "with", "yield", "in"
            ];
            var commonBuiltins = ["abs", "all", "any", "bin", "bool", "bytearray", "callable", "chr",
                "classmethod", "compile", "complex", "delattr", "dict", "dir", "divmod",
                "enumerate", "eval", "filter", "float", "format", "frozenset",
                "getattr", "globals", "hasattr", "hash", "help", "hex", "id",
                "int", "isinstance", "issubclass", "iter", "len",
                "list", "locals", "map", "max", "memoryview", "min", "next",
                "object", "oct", "open", "ord", "pow", "property", "range",
                "repr", "reversed", "round", "set", "setattr", "slice",
                "sorted", "staticmethod", "str", "sum", "super", "tuple",
                "type", "vars", "zip", "__import__", "NotImplemented",
                "Ellipsis", "__debug__",

                'group_by', 'filetype', 'paired_with', 'for_each', 'pattern', 'dynamic',
                'pattern', 'workdir', 'concurrent', 'docker_image', 'docker_file',
                'shared', 'skip', 'sigil', 'provides', 'input',
                'output'
            ];
            CodeMirror.registerHelper("hintWords", "sos", commonKeywords.concat(commonBuiltins));

            function top(state) {
                return state.scopes[state.scopes.length - 1];
            }

            CodeMirror.defineMode("sos", function(conf, parserConf) {
                var ERRORCLASS = "error";

                var singleDelimiters = parserConf.singleDelimiters || /^[\(\)\[\]\{\}@,:`=;\.]/;
                var doubleOperators = parserConf.doubleOperators || /^([!<>]==|<>|<<|>>|\/\/|\*\*)/;
                var doubleDelimiters = parserConf.doubleDelimiters || /^(\+=|\-=|\*=|%=|\/=|&=|\|=|\^=)/;
                var tripleDelimiters = parserConf.tripleDelimiters || /^(\/\/=|>>=|<<=|\*\*=)/;

                var hangingIndent = parserConf.hangingIndent || conf.indentUnit;

                var myKeywords = commonKeywords,
                    myBuiltins = commonBuiltins;
                if (parserConf.extra_keywords != undefined)
                    myKeywords = myKeywords.concat(parserConf.extra_keywords);

                if (parserConf.extra_builtins != undefined)
                    myBuiltins = myBuiltins.concat(parserConf.extra_builtins);

                var singleOperators = parserConf.singleOperators || /^[\+\-\*\/%\$&|\^~<>!@]/;
                var identifiers = parserConf.identifiers || /^[_A-Za-z\u00A1-\uFFFF][_A-Za-z0-9\u00A1-\uFFFF]*/;
                myKeywords = myKeywords.concat(["nonlocal", "False", "True", "None", "async", "await"]);
                myBuiltins = myBuiltins.concat(["ascii", "bytes", "exec", "print"]);
                var stringPrefixes = new RegExp("^(([rbuf]|(br))?('{3}|\"{3}|['\"]))", "i");
                var keywords = wordRegexp(myKeywords);
                var builtins = wordRegexp(myBuiltins);

                // tokenizers
                function tokenBase(stream, state) {
                    if (stream.sol()) state.indent = stream.indentation()
                    // Handle scope changes
                    if (stream.sol() && top(state).type == "py") {
                        var scopeOffset = top(state).offset;
                        if (stream.eatSpace()) {
                            var lineOffset = stream.indentation();
                            if (lineOffset > scopeOffset)
                                pushPyScope(state);
                            else if (lineOffset < scopeOffset && dedent(stream, state))
                                state.errorToken = true;
                            return null;
                        } else {
                            var style = tokenBaseInner(stream, state);
                            if (scopeOffset > 0 && dedent(stream, state))
                                style += " " + ERRORCLASS;
                            return style;
                        }
                    }
                    return tokenBaseInner(stream, state);
                }

                function tokenBaseInner(stream, state) {
                    if (stream.eatSpace()) return null;

                    var ch = stream.peek();

                    // Handle Comments
                    if (ch == "#") {
                        stream.skipToEnd();
                        return "comment";
                    }
                    // BO PENG
                    // handle shell command
                    if (state.beginningOfLine && stream.match(/![a-zA-Z]+/)) {
                        stream.next();
                        return "meta";
                    }
                    // handle magic
                    if (state.beginningOfLine && stream.match(/%[a-zA-Z]+/)) {
                        stream.next();
                        return "meta";
                    }

                    if (state.beginningOfLine && stream.match(/[a-zA-Z]+:/)) {
                        stream.next();
                        return "meta";
                    }
                    // handle section header
                    if (state.beginningOfLine && stream.match(/\[.*\]\s*$/)) {
                        stream.skipToEnd();
                        return "meta";
                    }
                    // Handle Number Literals
                    if (stream.match(/^[0-9\.]/, false)) {
                        var floatLiteral = false;
                        // Floats
                        if (stream.match(/^\d*\.\d+(e[\+\-]?\d+)?/i)) {
                            floatLiteral = true;
                        }
                        if (stream.match(/^\d+\.\d*/)) {
                            floatLiteral = true;
                        }
                        if (stream.match(/^\.\d+/)) {
                            floatLiteral = true;
                        }
                        if (floatLiteral) {
                            // Float literals may be "imaginary"
                            stream.eat(/J/i);
                            return "number";
                        }
                        // Integers
                        var intLiteral = false;
                        // Hex
                        if (stream.match(/^0x[0-9a-f]+/i)) intLiteral = true;
                        // Binary
                        if (stream.match(/^0b[01]+/i)) intLiteral = true;
                        // Octal
                        if (stream.match(/^0o[0-7]+/i)) intLiteral = true;
                        // Decimal
                        if (stream.match(/^[1-9]\d*(e[\+\-]?\d+)?/)) {
                            // Decimal literals may be "imaginary"
                            stream.eat(/J/i);
                            // TODO - Can you have imaginary longs?
                            intLiteral = true;
                        }
                        // Zero by itself with no other piece of number.
                        if (stream.match(/^0(?![\dx])/i)) intLiteral = true;
                        if (intLiteral) {
                            // Integer literals may be "long"
                            stream.eat(/L/i);
                            return "number";
                        }
                    }

                    // Handle Strings
                    if (stream.match(stringPrefixes)) {
                        state.tokenize = tokenStringFactory(stream.current());
                        return state.tokenize(stream, state);
                    }

                    // Handle operators and Delimiters
                    if (stream.match(tripleDelimiters) || stream.match(doubleDelimiters))
                        return "punctuation";

                    if (stream.match(doubleOperators) || stream.match(singleOperators))
                        return "operator";

                    if (stream.match(singleDelimiters))
                        return "punctuation";

                    if (state.lastToken == "." && stream.match(identifiers))
                        return "property";

                    if (stream.match(keywords) || stream.match(wordOperators))
                        return "keyword";

                    if (stream.match(builtins))
                        return "builtin";

                    if (stream.match(/^(self|cls)\b/))
                        return "variable-2";

                    if (stream.match(identifiers)) {
                        if (state.lastToken == "def" || state.lastToken == "class")
                            return "def";
                        return "variable";
                    }

                    // Handle non-detected items
                    stream.next();
                    return ERRORCLASS;
                }

                function tokenStringFactory(delimiter) {
                    while ("rub".indexOf(delimiter.charAt(0).toLowerCase()) >= 0)
                        delimiter = delimiter.substr(1);

                    var singleline = delimiter.length == 1;
                    var OUTCLASS = "string";

                    function tokenString(stream, state) {
                        while (!stream.eol()) {
                            stream.eatWhile(/[^'"\\]/);
                            if (stream.eat("\\")) {
                                stream.next();
                                if (singleline && stream.eol())
                                    return OUTCLASS;
                            } else if (stream.match(delimiter)) {
                                state.tokenize = tokenBase;
                                return OUTCLASS;
                            } else {
                                stream.eat(/['"]/);
                            }
                        }
                        if (singleline) {
                            if (parserConf.singleLineStringErrors)
                                return ERRORCLASS;
                            else
                                state.tokenize = tokenBase;
                        }
                        return OUTCLASS;
                    }
                    tokenString.isString = true;
                    return tokenString;
                }

                function pushPyScope(state) {
                    while (top(state).type != "py") state.scopes.pop()
                    state.scopes.push({
                        offset: top(state).offset + conf.indentUnit,
                        type: "py",
                        align: null
                    })
                }

                function pushBracketScope(stream, state, type) {
                    var align = stream.match(/^([\s\[\{\(]|#.*)*$/, false) ? null : stream.column() + 1
                    state.scopes.push({
                        offset: state.indent + hangingIndent,
                        type: type,
                        align: align
                    })
                }

                function dedent(stream, state) {
                    var indented = stream.indentation();
                    while (top(state).offset > indented) {
                        if (top(state).type != "sos") return true;
                        state.scopes.pop();
                    }
                    return top(state).offset != indented;
                }

                function tokenLexer(stream, state) {
                    if (stream.sol()) state.beginningOfLine = true;

                    var style = state.tokenize(stream, state);
                    var current = stream.current();

                    // Handle decorators
                    if (state.beginningOfLine && current == "@")
                        return stream.match(identifiers, false) ? "meta" : py3 ? "operator" : ERRORCLASS;



                    if (/\S/.test(current)) state.beginningOfLine = false;

                    //if ((style == "variable" || style == "builtin")
                    //    && state.lastToken == "meta")
                    //  style = "meta";

                    // Handle scope changes.
                    if (current == "pass" || current == "return")
                        state.dedent += 1;

                    if (current == "lambda") state.lambda = true;
                    if (current == ":" && !state.lambda && top(state).type == "py")
                        pushPyScope(state);

                    var delimiter_index = current.length == 1 ? "[({".indexOf(current) : -1;
                    if (delimiter_index != -1)
                        pushBracketScope(stream, state, "])}".slice(delimiter_index, delimiter_index + 1));

                    delimiter_index = "])}".indexOf(current);
                    if (delimiter_index != -1) {
                        if (top(state).type == current) state.indent = state.scopes.pop().offset - hangingIndent
                        else return ERRORCLASS;
                    }
                    if (state.dedent > 0 && stream.eol() && top(state).type == "py") {
                        if (state.scopes.length > 1) state.scopes.pop();
                        state.dedent -= 1;
                    }

                    return style;
                }

                var external = {
                    startState: function(basecolumn) {
                        return {
                            tokenize: tokenBase,
                            scopes: [{
                                offset: basecolumn || 0,
                                type: "sos",
                                align: null
                            }],
                            indent: basecolumn || 0,
                            lastToken: null,
                            lambda: false,
                            dedent: 0
                        };
                    },

                    token: function(stream, state) {
                        var addErr = state.errorToken;
                        if (addErr) state.errorToken = false;
                        var style = tokenLexer(stream, state);

                        if (style && style != "comment")
                            state.lastToken = (style == "keyword" || style == "punctuation") ? stream.current() : style;
                        if (style == "punctuation") style = null;

                        if (stream.eol() && state.lambda)
                            state.lambda = false;
                        return addErr ? style + " " + ERRORCLASS : style;
                    },

                    indent: function(state, textAfter) {
                        if (state.tokenize != tokenBase)
                            return state.tokenize.isString ? CodeMirror.Pass : 0;

                        var scope = top(state),
                            closing = scope.type == textAfter.charAt(0)
                        if (scope.align != null)
                            return scope.align - (closing ? 1 : 0)
                        else
                            return scope.offset - (closing ? hangingIndent : 0)
                    },

                    electricInput: /^\s*[\}\]\)]$/,
                    closeBrackets: {
                        triples: "'\""
                    },
                    lineComment: "#",
                    fold: "indent"
                };
                return external;
            });

            CodeMirror.defineMIME("text/x-sos", "sos");

        });

    }

    return {
        onload: onload
    }
})
