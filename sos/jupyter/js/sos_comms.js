// override the existing execute function by
        //
        // look for all input cells, find one that has prompt '*', which must
        // be the one that is being executed. Then, get the metadata and send
        // the kernel and cell index through the %softwith magic.
        //

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
                        // in addition, the softwidth command will send a "--list-kernel" request if
                        // the frontend is not correctly initialized, possibly because the kernel was
                        // not ready when the frontend sent the command `%listkernel`.
                        "%softwith " +
                        (window.kernel_updated ? "" : " --list-kernel ") +
                        " --default-kernel " + window.default_kernel +
                        " --cell-kernel " + cells[i].metadata.kernel +
                        " --cell " + i.toString() + "\n" + code,
                        callbacks, options)
                }
            }
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
                        // the real kernel name back to kernel (%softwith and metadata).
                        //
                        // there are two kinds of messages from my_execute
                        // 1. cell_idx: kernel
                        //     the kernel used for the cell with source
                        // 2. None: kernel
                        //     the kernel for the new cell

                        var data = msg.content.data;
                        // console.log(data);

                        if (data[0] instanceof Array) {
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
                        } else {
                            // update the cells when the notebook is being opened.
                            if (data[0] == null) {
                                // we also set a global kernel to be used for new cells
                                $('#kernel_selector').val(DisplayName[data[1]]);
                                // a side effect of change is cells without metadata kernel info will change background
                                $('#kernel_selector').change();
                            } else {
                                // get cell from passed cell index, which was sent through the
                                // %softwith magic
                                var cell = IPython.notebook.get_cell(data[0]);
                                if (cell.metadata.kernel != KernelName[data[1]]) {
                                    cell.metadata.kernel = KernelName[data[1]];
                                    // set meta information
                                    changeStyleOnKernel(cell, data[1])
                                }
                            }
                        }
                    });
                }
            );
        }

        

        function wrap_execute() {
                if (!window.kernel_updated)
                    IPython.notebook.kernel.execute('%softwith --list-kernel',
                        [], {'silent': true, 'store_history': false});
                // override kernel execute with the wrapper.
                IPython.notebook.kernel.orig_execute = IPython.notebook.kernel.execute
                IPython.notebook.kernel.execute = my_execute
            }