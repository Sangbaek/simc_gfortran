import subprocess

def ndelta():
  with open('scripts/ndelta_kinema_output.txt', 'r') as kinema_file:
    kinema_lines = kinema_file.readlines()

  with open('infiles/ndelta_template.inp', 'r') as template_file:
    template_content = template_file.read()

  runstrings = '#!/bin/tcsh\n\ngo_ndelta_vcs\n'

  for kinema_setting in kinema_lines:
    kin_index          = kinema_setting.split()[0]
    custom_beam_charge = kinema_setting.split()[1]
    custom_ebeam       = kinema_setting.split()[2]
    custom_ep          = kinema_setting.split()[3]
    custom_etheta      = kinema_setting.split()[4]
    custom_pp          = kinema_setting.split()[5]
    custom_ptheta      = kinema_setting.split()[6]
    print(kin_index)

    customized_content = template_content.replace('custom_beam_charge', str(custom_beam_charge))
    customized_content = customized_content.replace('custom_ebeam', str(custom_ebeam))
    customized_content = customized_content.replace('custom_ep', str(custom_ep))
    customized_content = customized_content.replace('custom_etheta', str(custom_etheta))
    customized_content = customized_content.replace('custom_pp', str(custom_pp))
    customized_content = customized_content.replace('custom_ptheta', str(custom_ptheta))


    for custom_random_seed in range(1000, 1001):
      customized_content_with_random_seed = customized_content.replace('custom_random_seed', str(custom_random_seed))
      customized_filename = 'infiles/nd26_{}_e{}.inp'.format(kin_index, custom_ebeam)
      with open(customized_filename, 'w') as custom_file:
        custom_file.write(customized_content_with_random_seed)
      runstrings = runstrings + 'echo {} {}\n'.format(kin_index, customized_filename)
      runstrings = runstrings + './run_simc nd26_{}_e{}\n'.format(kin_index, custom_ebeam)
      runstrings = runstrings + 'h2root worksim/nd26_{}_e{}.rzdat\n'.format(kin_index, custom_ebeam)
  with open('/work/hallc/alphaE/sangbaek/simc_gfortran/scripts/run_infiles_ndelta.csh', 'w') as runfile:
    runfile.write(runstrings)

  subprocess.run(['chmod', '+x', '/work/hallc/alphaE/sangbaek/simc_gfortran/scripts/run_infiles_ndelta.csh'])

def vcs():
  with open('scripts/vcs_kinema_output.txt', 'r') as kinema_file:
    kinema_lines = kinema_file.readlines()

  with open('infiles/vcs_template.inp', 'r') as template_file:
    template_content = template_file.read()

  runstrings = '#!/bin/tcsh\n\ngo_ndelta_vcs\n'

  for kinema_setting in kinema_lines:
    kin_index          = kinema_setting.split()[0]
    custom_beam_charge = kinema_setting.split()[1]
    custom_ebeam       = kinema_setting.split()[2]
    custom_ep          = kinema_setting.split()[3]
    custom_etheta      = kinema_setting.split()[4]
    custom_pp          = kinema_setting.split()[5]
    custom_ptheta      = kinema_setting.split()[6]
    print(kin_index)

    customized_content = template_content.replace('custom_beam_charge', str(custom_beam_charge))
    customized_content = customized_content.replace('custom_ebeam', str(custom_ebeam))
    customized_content = customized_content.replace('custom_ep', str(custom_ep))
    customized_content = customized_content.replace('custom_etheta', str(custom_etheta))
    customized_content = customized_content.replace('custom_pp', str(custom_pp))
    customized_content = customized_content.replace('custom_ptheta', str(custom_ptheta))
    customized_content = customized_content.replace('custom_kin', 'v2{}'.format(kin_index[:3]))


    for custom_random_seed in range(1000, 1001):
      customized_content_with_random_seed = customized_content.replace('custom_random_seed', str(custom_random_seed))
      customized_filename = 'infiles/vcs2_{}_e{}.inp'.format(kin_index, custom_ebeam)
      with open(customized_filename, 'w') as custom_file:
        custom_file.write(customized_content_with_random_seed)
      runstrings = runstrings + 'echo {} {}\n'.format(kin_index, customized_filename)
      runstrings = runstrings + './run_simc vcs2_{}_e{}\n'.format(kin_index, custom_ebeam)
      runstrings = runstrings + 'h2root worksim/vcs2_{}_e{}.rzdat\n'.format(kin_index, custom_ebeam)
  with open('/work/hallc/alphaE/sangbaek/simc_gfortran/scripts/run_infiles_vcs.csh', 'w') as runfile:
    runfile.write(runstrings)

  subprocess.run(['chmod', '+x', '/work/hallc/alphaE/sangbaek/simc_gfortran/scripts/run_infiles_vcs.csh'])

def main():
  ndelta()
  vcs()

if __name__ == '__main__':
  

  main()
