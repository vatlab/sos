# sample ipython_config.py for SoS environment
c = get_config()

c.TerminalIPythonApp.display_banner = True
c.InteractiveShellApp.log_level = 20
c.InteractiveShellApp.extensions = [
    'sos_magic'
]
c.InteractiveShell.autoindent = True
c.InteractiveShell.colors = 'LightBG'
c.InteractiveShell.confirm_exit = False
c.InteractiveShell.deep_reload = True

c.PromptManager.in_template  = 'In [\#]: '
c.PromptManager.in2_template = '   .\D.: '
c.PromptManager.out_template = 'Out[\#]: '
c.PromptManager.justify = True

c.PrefilterManager.multi_line_specials = True

c.AliasManager.user_aliases = [
]
