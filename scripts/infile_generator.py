import subprocess

def main():
  with open('scripts/kinema_output.txt', 'r') as kinema_file:
    kinema_lines = kinema_file.readlines()

  with open('infiles/ndelta_template.inp', 'r') as template_file:
    template_content = template_file.read()

  runstrings = '#!/bin/tcsh\n\ncd /volatile/hallc/alphaE/sangbaek/simc_practice/simc_gfortran/\nmodule use /scigroup/cvmfs/scicomp/sw/el9/modulefiles\nmodule load root/6.30.06-gcc11.4.0\n'

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


    for custom_random_seed in range(1411, 1420):
      customized_content_with_random_seed = customized_content.replace('custom_random_seed', str(custom_random_seed))
      customized_filename = 'infiles/nd26_{}_e{}_rs{}.inp'.format(kin_index, custom_ebeam, custom_random_seed)
      with open(customized_filename, 'w') as custom_file:
        custom_file.write(customized_content_with_random_seed)
      runstrings = runstrings + 'echo {} {}\n'.format(kin_index, customized_filename)
      runstrings = runstrings + './run_simc nd26_{}_e{}_rs{}\n'.format(kin_index, custom_ebeam, custom_random_seed)
      runstrings = runstrings + 'h2root worksim/nd26_{}_e{}_rs{}.rzdat\n'.format(kin_index, custom_ebeam, custom_random_seed)
  with open('run_infiles.csh', 'w') as runfile:
    runfile.write(runstrings)

  subprocess.run(['chmod', '+x', 'run_infiles.csh'])

if __name__ == '__main__':
  

  main()
