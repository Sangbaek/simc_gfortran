with open('kinema_output.txt', 'r') as kinema_file:
  kinema_lines = kinema_file.readlines()

with open('infiles/ndelta_template.inp', 'r') as template_file:
  template_content = template_file.read()


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
  customimzed_content = customized_content.replace('custom_ep', str(custom_ep))
  customized_content = customized_content.replace('custom_etheta', str(custom_etheta))
  customized_content = customized_content.replace('custom_pp', str(custom_pp))
  customized_content = customized_content.replace('custom_ptheta', str(custom_ptheta))


  for custom_random_seed in range(1410, 1420):
    customized_content = customized_content.replace('custom_random_seed', str(custom_random_seed))
    with open('infiles/ndelta_2026_kin{}_ebeam_{}_randomseed_{}.inp'.format(kin_index, custom_ebeam, custom_random_seed), 'w') as custom_file:
      custom_file.write(customized_content)

