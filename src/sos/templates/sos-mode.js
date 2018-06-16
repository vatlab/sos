/**
 * Copyright (c) Bo Peng and UT MD Anderson Cancer Center
 * Distributed under the terms of the Modified BSD License.
 **/

(function(mod) {
  if (typeof exports == "object" && typeof module == "object") // CommonJS
    mod(require("codemirror/lib/codemirror"),
      require("codemirror/mode/python/python"),
      require("codemirror/mode/r/r"),
      require("codemirror/mode/octave/octave"),
      require("codemirror/mode/ruby/ruby"),
      require("codemirror/mode/sas/sas"),
      require("codemirror/mode/javascript/javascript"),
      require("codemirror/mode/shell/shell"),
      require("codemirror/mode/julia/julia"),
      require("codemirror/mode/markdown/markdown"));
  else if (typeof define == "function" && define.amd) // AMD
    define(["codemirror/lib/codemirror", "codemirror/mode/python/python",
      "codemirror/mode/markdown/markdown", "codemirror/mode/r/r"
    ], mod);
  else // Plain browser env
    mod(CodeMirror);
})(function(CodeMirror) {
  "use strict";

  var sosKeywords = ["input", "output", "depends", "parameter"];
  var sosActionWords = ["script", "download", "run", "bash", "sh", "csh",
    "tcsh", "zsh", "python", "python2", "python3", "R", "node", "julia",
    "matlab", "octave", "ruby", "perl", "report", "pandoc", "docker_build",
    "Rmarkdown"
  ];
  var sosMagicWords = ['cd', 'capture', 'clear', 'debug', 'dict', 'expand', 'get',
    'matplotlib', 'paste', 'preview', 'pull', 'push', 'put', 'render',
    'rerun', 'run', 'save', 'sandbox', 'set', 'sessioninfo', 'sosrun',
    'sossave', 'shutdown', 'taskinfo', 'tasks', 'toc', 'use', 'with'
  ]
  var sosFunctionWords = ["sos_run", "logger", "get_output"];

  var hintWords = sosKeywords.concat(sosActionWords).concat(sosFunctionWords)
    .concat(sosMagicWords);

  var sosDirectives = sosKeywords.map(x => x + ":");
  var sosActions = sosActionWords.map(x => x + ":");
  var sosMagics = sosMagicWords.map(x => '%' + x);

  // hint word for SoS mode
  CodeMirror.registerHelper("hintWords", "sos", hintWords);

  var modeMap = {
    'sos': null,
    'python': {
      name: 'python',
      version: 3
    },
    'python2': {
      name: 'python',
      version: 2
    },
    'python3': {
      name: 'python',
      version: 3
    },
    'r': 'r',
    'report': 'markdown',
    'pandoc': 'markdown',
    'download': 'markdown',
    'markdown': 'markdown',
    'ruby': 'ruby',
    'sas': 'sas',
    'bash': 'shell',
    'sh': 'shell',
    'julia': 'julia',
    'run': 'shell',
    'javascript': 'javascript',
    'typescript': {
      name: "javascript",
      typescript: true
    },
    'octave': 'octave',
    'matlab': 'octave',
  }

  function findMode(mode) {
    if (mode in modeMap) {
      return modeMap[mode];
    }
    return null;
  }

  function markExpr(python_mode) {
    return {
      startState: function() {
        return {
          in_python: false,
          sigil: null,
          matched: true,
          python_state: CodeMirror.startState(python_mode),
        };
      },

      copyState: function(state) {
        return {
          in_python: state.in_python,
          sigil: state.sigil,
          matched: state.matched,
          python_state: CodeMirror.copyState(python_mode, state.python_state)
        };
      },

      token: function(stream, state) {
        if (state.in_python) {
          if (stream.match(state.sigil.right)) {
            state.in_python = false;
            state.python_state = CodeMirror.startState(python_mode);
            return "sos-sigil";
          }
          let it = null;
          try {
            it = python_mode.token(stream, state.python_state);
          } catch (error) {
            return "sos-interpolated error" + (state.matched ? "" : " sos-unmatched");
          }
          if (it == 'variable' || it == 'builtin') {
            let ct = stream.current();
            // warn users in the use of input and output in {}
            if (ct === 'input' || ct === 'output')
              it += ' error';
          }
          return (it ? ("sos-interpolated " + it) : "sos-interpolated") + (state.matched ? "" : " sos-unmatched");
        } else {
          // remove the double brace case, the syntax highlighter
          // does not have to worry (highlight) }}, although it would
          // probably mark an error for single }
          if (state.sigil.left === '{' && stream.match(/\{\{/))
            return null;
          if (stream.match(state.sigil.left)) {
            state.in_python = true;
            // let us see if there is any right sigil till the end of the editor.
            try {
              let rest = stream.string.slice(stream.pos);
              if (!rest.includes(state.sigil.right)) {
                state.matched = false;
                for (let idx = 1; idx < 5; ++idx) {
                  if (stream.lookAhead(idx).includes(state.sigil.right)) {
                    state.matched = true;
                    break;
                  }
                }
              }
            } catch (error) {
              // only codemirror 5.27.0 supports this function
            }
            return "sos-sigil" + (state.matched ? "" : " sos-unmatched");
          }
          while (stream.next() && !stream.match(state.sigil.left, false)) {}
          return null;
        }
      }
    }
  }

  CodeMirror.defineMode("sos", function(conf, parserConf) {
    let sosPythonConf = {};
    for (let prop in parserConf) {
      if (parserConf.hasOwnProperty(prop)) {
        sosPythonConf[prop] = parserConf[prop];
      }
    }
    sosPythonConf.name = 'python';
    sosPythonConf.version = 3;
    sosPythonConf.extra_keywords = sosActionWords.concat(sosFunctionWords);
    // this is the SoS flavored python mode with more identifiers
    var base_mode = null;
    if ('base_mode' in parserConf && parserConf.base_mode && parserConf.base_mode !== 'sos' && parserConf.base_mode !== 'SoS') {
      let mode = findMode(parserConf.base_mode.toLowerCase());
      if (mode) {
        base_mode = CodeMirror.getMode(conf, mode);
      } else {
        console.log(`No base mode is found for ${parserConf.base_mode}. Python mode used.`);
      }
    } else if ('base_mode' in conf && conf.base_mode && conf.base_mode !== 'sos' && conf.base_mode !== 'SoS') {
      let mode = findMode(conf.base_mode.toLowerCase());
      if (mode) {
        base_mode = CodeMirror.getMode(conf, mode);
      } else {
        console.log(`No base mode is found for ${conf.base_mode}. Python mode used.`);
      }
    }
    // if there is a user specified base mode, this is the single cell mode
    if (base_mode) {
      var python_mode = CodeMirror.getMode({}, {
        name: 'python',
        version: 3
      });
      var overlay_mode = markExpr(python_mode);
      return {
        startState: function() {
          return {
            sos_mode: true,
            base_state: CodeMirror.startState(base_mode),
            overlay_state: CodeMirror.startState(overlay_mode),
            // for overlay
            basePos: 0,
            baseCur: null,
            overlayPos: 0,
            overlayCur: null,
            streamSeen: null
          };
        },

        copyState: function(state) {
          return {
            sos_mode: state.sos_mode,
            base_state: CodeMirror.copyState(base_mode, state.base_state),
            overlay_state: CodeMirror.copyState(overlay_mode, state.overlay_state),
            // for overlay
            basePos: state.basePos,
            baseCur: null,
            overlayPos: state.overlayPos,
            overlayCur: null
          };
        },

        token: function(stream, state) {
          if (state.sos_mode) {
            if (stream.sol()) {
              let sl = stream.peek();
              if (sl == '!') {
                stream.skipToEnd();
                return "meta";
              } else if (sl == '#') {
                stream.skipToEnd();
                return 'comment'
              }
              for (var i = 0; i < sosMagics.length; i++) {
                if (stream.match(sosMagics[i])) {
                  if (sosMagics[i] === "%expand") {
                    if (stream.eol() || stream.match(/\s*$/, false)) {
                      state.overlay_state.sigil = {
                        'left': '{',
                        'right': '}'
                      }
                    } else {
                      let found = stream.match(/\s+(\S+)\s+(\S+)$/, false);
                      if (found) {
                        state.overlay_state.sigil = {
                          'left': found[1],
                          'right': found[2]
                        }
                      } else {
                        state.overlay_state.sigil = null;
                      }
                    }
                  }
                  // the rest of the lines will be processed as Python code
                  return "meta";
                }
              }
              state.sos_mode = false;
            } else {
              stream.skipToEnd();
              return null;
            }
          }

          if (state.overlay_state.sigil) {
            if (stream != state.streamSeen ||
              Math.min(state.basePos, state.overlayPos) < stream.start) {
              state.streamSeen = stream;
              state.basePos = state.overlayPos = stream.start;
            }

            if (stream.start == state.basePos) {
              state.baseCur = base_mode.token(stream, state.base_state);
              state.basePos = stream.pos;
            }
            if (stream.start == state.overlayPos) {
              stream.pos = stream.start;
              state.overlayCur = overlay_mode.token(stream, state.overlay_state);
              state.overlayPos = stream.pos;
            }
            stream.pos = Math.min(state.basePos, state.overlayPos);

            // state.overlay.combineTokens always takes precedence over combine,
            // unless set to null
            return state.overlayCur ? state.overlayCur : state.baseCur;
          } else {
            return base_mode.token(stream, state.base_state);
          }
        },

        indent: function(state, textAfter) {
          // inner indent
          if (!state.sos_mode) {
            if (!base_mode.indent) return CodeMirror.Pass;
            return base_mode.indent(state.base_state, textAfter);
          } else {
            // sos mode has no indent
            return 0;
          }
        },

        innerMode: function(state) {
          return state.sos_mode ? {
            state: state.base_state,
            mode: base_mode
          } : null;
        },

        lineComment: "#",
        fold: "indent"
      };
    } else {
      // this is SoS mode
      base_mode = CodeMirror.getMode(conf, sosPythonConf);
      overlay_mode = markExpr(base_mode);
      return {
        startState: function() {
          return {
            sos_state: null,
            base_state: CodeMirror.startState(base_mode),
            overlay_state: CodeMirror.startState(overlay_mode),
            inner_mode: null,
            inner_state: null,
            // for overlay
            basePos: 0,
            baseCur: null,
            overlayPos: 0,
            overlayCur: null,
            streamSeen: null
          };
        },

        copyState: function(state) {
          return {
            sos_state: state.sos_state,
            base_state: CodeMirror.copyState(base_mode, state.base_state),
            overlay_state: CodeMirror.copyState(overlay_mode, state.overlay_state),
            inner_mode: state.inner_mode,
            inner_state: state.inner_mode && CodeMirror.copyState(state.inner_mode, state.inner_state),
            // for overlay
            basePos: state.basePos,
            baseCur: null,
            overlayPos: state.overlayPos,
            overlayCur: null
          };
        },

        token: function(stream, state) {
          if (stream.sol()) {
            let sl = stream.peek();
            if (sl == '[') {
              // header, move to the end
              if (stream.match(/^\[.*\]$/, false)) {
                // if there is no :, the easy case
                if (stream.match(/^\[[^:]*\]$/)) {
                  // reset state
                  state.sos_state = null;
                  state.inner_mode = null;
                  return "header line-section-header";
                } else {
                  // match up to :
                  stream.match(/^\[[^:]*:/);
                  state.sos_state = 'header_option';
                  return "header line-section-header";
                }
              }
            } else if (sl == '!') {
              stream.eatWhile(/\S/);
              return "meta";
            } else if (sl == '#') {
              stream.skipToEnd();
              return "comment";
            } else if (sl == '%') {
              stream.eatWhile(/\S/);
              return "meta";
            } else if (state.sos_state && state.sos_state.startsWith('entering')) {
              // the second parameter is starting column
              let mode = findMode(state.sos_state.slice(9).toLowerCase());
              state.inner_mode = CodeMirror.getMode(conf, mode);
              state.inner_state = CodeMirror.startState(state.inner_mode, stream.indentation());
              state.sos_state = null;
            }
            for (var i = 0; i < sosDirectives.length; i++) {
              if (stream.match(sosDirectives[i])) {
                // the rest of the lines will be processed as Python code
                state.sos_state = 'directive_option'
                return "keyword strong";
              }
            }
            for (var i = 0; i < sosActions.length; i++) {
              if (stream.match(sosActions[i])) {
                // switch to submode?
                if (stream.eol()) {
                  // really
                  let mode = findMode(stream.current().slice(0, -1).toLowerCase());
                  if (mode) {
                    state.sos_state = "entering " + stream.current().slice(0, -1);
                  } else {
                    state.sos_state = 'unknown_language';
                  }
                } else {
                  state.sos_state = 'start ' + stream.current().slice(0, -1);
                }
                state.overlay_state.sigil = null;
                return "builtin strong";
              }
            }
            // if unknown action
            if (stream.match(/\w+:/)) {
              state.overlay_state.sigil = null;
              state.sos_state = 'start ' + stream.current().slice(0, -1);
              return "builtin strong";
            }
          } else if (state.sos_state == 'header_option') {
            // stuff after :
            if (stream.peek() == ']') {
              // move next
              stream.next();
              // ] is the last char
              if (stream.eol()) {
                state.sos_state = null;
                state.inner_mode = null;
                return "header line-section-header";
              } else {
                stream.backUp(1);
                let it = base_mode.token(stream, state.base_state);
                return it ? it + ' sos-option' : null;
              }
            } else {
              let it = base_mode.token(stream, state.base_state);
              return it ? it + ' sos-option' : null;
            }
          } else if (state.sos_state == 'directive_option') {
            // stuff after input:, R: etc
            if (stream.peek() == ',') {
              // move next
              stream.next();
              // , is the last char, continue option line
              if (stream.eol()) {
                stream.backUp(1);
                let it = base_mode.token(stream, state.base_state);
                return it ? it + ' sos-option' : null;
              }
              stream.backUp(1);
            } else if (stream.eol()) {
              // end of line stops option mode
              state.sos_state = null;
              state.inner_mode = null;
            }
            let it = base_mode.token(stream, state.base_state);
            return it ? it + ' sos-option' : null;
          } else if (state.sos_state && state.sos_state.startsWith("start ")) {
            // try to understand option expand=
            if (stream.match(/expand\s*=\s*True/, false)) {
              // highlight {}
              state.overlay_state.sigil = {
                'left': '{',
                'right': '}'
              }
            } else {
              let found = stream.match(/expand\s*=\s*"(\S+) (\S+)"/, false);
              if (!found)
                found = stream.match(/expand\s*=\s*'(\S+) (\S+)'/, false);
              if (found) {
                state.overlay_state.sigil = {
                  'left': found[1],
                  'right': found[2]
                }
              }
            }
            let token = base_mode.token(stream, state.base_state);
            // if it is end of line, ending the starting switch mode
            if (stream.eol() && stream.peek() !== ',') {
              // really
              let mode = findMode(state.sos_state.slice(6).toLowerCase());
              if (mode) {
                state.sos_state = "entering " + state.sos_state.slice(6);
              } else {
                state.sos_state = 'unknown_language';
              }
            }
            return token + ' sos-option';
          }
          // can be start of line but not special
          if (state.sos_state == 'unknown_language') {
            // we still handle {} in no man unknown_language
            if (state.overlay_state.sigil) {
              return overlay_mode.token(stream, state.overlay_state);
            } else {
              stream.skipToEnd();
              return null;
            }
          } else if (state.inner_mode) {
            let it = 'sos_script ';
            if (!state.overlay_state.sigil) {
              let st = state.inner_mode.token(stream, state.inner_state);
              return st ? it + st : null;
            } else {
              // overlay mode, more complicated
              if (stream != state.streamSeen ||
                Math.min(state.basePos, state.overlayPos) < stream.start) {
                state.streamSeen = stream;
                state.basePos = state.overlayPos = stream.start;
              }

              if (stream.start == state.basePos) {
                state.baseCur = state.inner_mode.token(stream, state.inner_state);
                state.basePos = stream.pos;
              }
              if (stream.start == state.overlayPos) {
                stream.pos = stream.start;
                state.overlayCur = overlay_mode.token(stream, state.overlay_state);
                state.overlayPos = stream.pos;
              }
              stream.pos = Math.min(state.basePos, state.overlayPos);
              // state.overlay.combineTokens always takes precedence over combine,
              // unless set to null
              return (state.overlayCur ? state.overlayCur : state.baseCur) + " sos-script";
            }
          } else {
            return base_mode.token(stream, state.base_state);
          }
        },

        indent: function(state, textAfter) {
          // inner indent
          if (state.inner_mode) {
            if (!state.inner_mode.indent) return CodeMirror.Pass;
            return state.inner_mode.indent(state.inner_mode, textAfter) + 2;
          } else {
            return base_mode.indent(state.base_state, textAfter);
          }
        },

        innerMode: function(state) {
          return state.inner_mode ? null : {
            state: state.base_state,
            mode: base_mode
          };
        },

        lineComment: "#",
        fold: "indent",
        electricInput: /^\s*[\}\]\)]$/,
      };
    };
  }, "python");

  CodeMirror.defineMIME("text/x-sos", "sos");
});
