var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SPM","category":"page"},{"location":"#SPM","page":"Home","title":"SPM","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SPM, a package for Statistical Process Monitoring.","category":"page"},{"location":"#Package-features","page":"Home","title":"Package features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Standard control charts\nUnivariate Shewhart, EWMA, AEWMA, and CUSUM control charts;\nMEWMA and MAEWMA control charts;\nControl limits\nClassical one-sided and two-sided control limits;\nSupport for multi-chart combinations;\nDynamic control limits based on bootstrap and permutation methods;\nOptimization methods\nState-of-the-art methods for estimating control limits;\nOptimization of control chart parameters against user-defined out-of-control scenarios;\nExtensibility to user-made control statistics\nUsers only need to define the behaviour of the control statistic (struct and sequential update function), everything else is taken care of by the package.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SPM]","category":"page"},{"location":"#SPM.get_limit-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_limit","text":"get_limit(CH::AbstractChart)\n\nGet the control limit of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.get_limit_value-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_limit_value","text":"get_limit_value(CH::AbstractChart)\n\nGet the control limit value of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.get_maxrl-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_maxrl","text":"get_maxrl(CH::AbstractChart)\n\nGet the maximum run length of the control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.get_nominal-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_nominal","text":"get_nominal(CH::AbstractChart)\nget_nominal_value(CH::AbstractChart)\n\nGet the nominal properties of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.get_param-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_param","text":"get_param(CH::AbstractChart)\nset_param!(CH::AbstractChart, par)\n\nGet and set the parameters of the control chart statistic.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.get_phase1-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_phase1","text":"get_phase1(CH::AbstractChart)\n\nGet the Phase 1 information of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.get_statistic-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_statistic","text":"get_statistic(CH::AbstractChart)\n\nGet the statistic of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.get_value-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.get_value","text":"get_value(CH::AbstractChart)\n\nGet the current value of the control chart statistic.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.is_IC-Tuple{SPM.AbstractChart}","page":"Home","title":"SPM.is_IC","text":"is_IC(CH::AbstractChart)\nis_OC(CH::AbstractChart)\n\nCheck whether the control chart is in control or out of control.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.new_data-Union{Tuple{Phase1Data{T}}, Tuple{T}} where T<:(AbstractVector)","page":"Home","title":"SPM.new_data","text":"new_data(P1::Phase1Data)\n\nGenerates a new observation via bootstrap, based on the observed Phase I data.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.set_limit!-Union{Tuple{C}, Tuple{LIM}, Tuple{C, LIM}} where {LIM<:SPM.AbstractLimit, C<:SPM.AbstractChart}","page":"Home","title":"SPM.set_limit!","text":"function set_limit!(CH::AbstractChart, limit::AbstractLimit)\n\nSet the control limit of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.set_nominal!-Union{Tuple{C}, Tuple{N}, Tuple{C, N}} where {N<:SPM.NominalProperties, C<:SPM.AbstractChart}","page":"Home","title":"SPM.set_nominal!","text":"set_nominal!(CH::AbstractChart, nominal::NominalProperties)\n\nSet the nominal properties of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.set_parameter!-Union{Tuple{C}, Tuple{C, Any}} where C<:SPM.AbstractChart","page":"Home","title":"SPM.set_parameter!","text":"function set_parameter!(CH::C, par)\n\nSet the parameters of the control chart statistic.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.set_phase1!-Union{Tuple{C}, Tuple{PH1}, Tuple{C, PH1}} where {PH1<:SPM.AbstractPhase1, C<:SPM.AbstractChart}","page":"Home","title":"SPM.set_phase1!","text":"set_phase1!(CH::AbstractChart, phase1::AbstractPhase1)\n\nSet the Phase 1 information of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.set_statistic!-Union{Tuple{C}, Tuple{STAT}, Tuple{C, STAT}} where {STAT<:SPM.AbstractStatistic, C<:SPM.AbstractChart}","page":"Home","title":"SPM.set_statistic!","text":"function set_statistic!(CH::AbstractChart, statistic::AbstractStatistic)\n\nSet the statistic of a control chart.\n\n\n\n\n\n","category":"method"},{"location":"#SPM.update_chart!-Tuple{SPM.AbstractChart, Any}","page":"Home","title":"SPM.update_chart!","text":"update_chart!(CH::AbstractChart, x)\n\nUpdate the control chart using a new observation x.\n\n\n\n\n\n","category":"method"}]
}
