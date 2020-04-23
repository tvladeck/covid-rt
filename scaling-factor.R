ACTUAL_ONSETS = 100
predicted_percentage = .045
actual_percentage = .06

accounted_for_onsets = actual_percentage * ACTUAL_ONSETS
we_scale_it_up_to = accounted_for_onsets / predicted_percentage
our_model_predicted = ACTUAL_ONSETS * predicted_percentage

# likelihood of observation according to model with raw counts
dpois(round(we_scale_it_up_to), ACTUAL_ONSETS)

# likelihood of observation according to model with scaled counts
dpois(accounted_for_onsets, our_model_predicted)
