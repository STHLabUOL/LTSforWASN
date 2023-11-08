function [K_Nominator,K_Denominator] = Controller_Choice(ControllerStructure, Tf, B_NumFr)
%Controller_Choice returns the Nominator and Denominator of the
%IMC-Controller for given parameter
%   Controller_Choice returns the Nominator and Denominator of the
%   IMC-Controller for given parameter
%   Parameter:
%   ControllerStructure: 'exact_IMC', 'approx_IMC', 'PIT1','PI'
%   Tf: '1','3','8'
 %Begin controller choice
    switch(ControllerStructure)
        case 'exact_IMC'
            switch B_NumFr
                case 7
                    switch Tf                     
                        case '1'
                            %Tf=1
                            K_Nominator = [0 0 21.4992496852606 -96.9616160805257 173.058210341506 -152.396420929236 65.9684267297528 -11.1678497467579 0 0 0 0];
                            K_Denominator = [1 -3.74970675828926 5.25152320728487 -3.25392613970197 0.752109690706334 1.00460604502252e-14 1.01859616086276e-14 -1.00000000000002 3.74970675828927 -5.25152320728486 3.25392613970196 -0.752109690706347];
                        case '3'
                            %Tf=3
                            K_Nominator = [0 0 2.59844540902351 -13.0961648614785 27.1272503811235 -29.5045369044942 17.7351333058324 -5.57550556942516 0.715378239418345 0 0 0 0];
                            K_Denominator = [1 -4.43646149417387 7.79064445949201 -6.75277178245544 2.87945616313032 -0.480867345993017 -2.31639508323559e-14 -0.999999999999970 4.43646149417386 -7.79064445949203 6.75277178245547 -2.87945616313034 0.480867345993030];
                        case '8'
                            %Tf=8
                            K_Nominator = [0 0 0.375232719232229 -1.89117290493043 3.91735454224061 -4.26065045427354 2.56107065908860 -0.805139915060232 0.103305353702769 0 0 0 0];
                            K_Denominator = [1 -4.48825464011057 7.97295633318921 -6.98919195571274 3.01253347230016 -0.508043209666078 2.29986166004310e-14 -0.999999999999974 4.48825464011051 -7.97295633318913 6.98919195571266 -3.01253347230008 0.508043209666013];
                        otherwise
                        error('No appropriate filter time constant Tf ={1,3,8} (parDXCPPhaT_cl.Tf)')
                    end
                case 39
                    switch Tf                           
                        case '1'
                            %Tf=1
                            K_Nominator = [0 0 119.781533960738 -603.698931162118 1250.49525824331 -1360.08194628558 817.543237581085 -257.016217232503 32.9770648950688 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
                            K_Denominator = [1 -4.27970675828930 7.23886778917839 -6.03723343956366 2.47669054475069 -0.398618136080176 7.85588887548677e-12 5.08885347370001e-12 -3.97289657503666e-11 -8.41175712890986e-12 6.34899122884490e-10 -3.25882329943958e-09 1.10495700245870e-08 -2.98328455933113e-08 6.79748603070510e-08 -1.33065590036055e-07 2.24147993874308e-07 -3.22964797483388e-07 3.93160308832986e-07 -3.94801463667217e-07 3.09386454105976e-07 -1.57102960636466e-07 -1.04017896259129e-08 1.38860624275881e-07 -1.99588375355763e-07 1.96747643166095e-07 -1.54165708078269e-07 9.79660845420690e-08 -4.72518092894676e-08 1.22529647801822e-08 5.08315750028491e-09 -9.16784322548733e-09 6.87883004420400e-09 -3.54081303073527e-09 1.37108897011315e-09 -4.23572203262735e-10 1.18080187641003e-10 -3.67098051878335e-11 1.34064782615642e-11 -1.00000000000474 4.27970675829069 -7.23886778917866 6.03723343956322 -2.47669054474853 0.398618136074394];
                        case '3'
                            %Tf=3
                            K_Nominator = [0 0 14.4770529931310 -72.9643470853799 151.137537837688 -164.382419896468 98.8100284182094 -31.0635310296545 3.98567876247364 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
                            K_Denominator = [1 -4.43646149417388 7.79064445949216 -6.75277178245642 2.87945616313648 -0.480867346027663 1.61961889494995e-10 -6.20813098799483e-10 1.95857992299508e-09 -5.12249162404791e-09 1.11672765689033e-08 -2.02645564699932e-08 3.00689026105831e-08 -3.42865815267950e-08 2.31418616815241e-08 1.25923137461141e-08 -7.43596551487134e-08 1.49270355862997e-07 -2.09948179956538e-07 2.23750957165103e-07 -1.68536845216049e-07 4.72945219427771e-08 1.07681256656998e-07 -2.47086287912151e-07 3.29295473444532e-07 -3.41621880016813e-07 3.02085475372409e-07 -2.41567784213751e-07 1.83102549963585e-07 -1.34435430808630e-07 9.42644257325192e-08 -6.08947342086718e-08 3.49782121257575e-08 -1.74137658325425e-08 7.37856211653990e-09 -2.61301513831198e-09 7.51208310312761e-10 -1.64632325957237e-10 2.26858754656628e-11 -0.999999999999750 4.43646149417274 -7.79064445949171 6.75277178245543 -2.87945616313035 0.480867345993031];
                        case '8'
                            %Tf=8
                            K_Nominator = [0 0 2.09058229286527 -10.5365347560410 21.8252610210549 -23.7379096738098 14.2688222434937 -4.48577952676418 0.575558399201145 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
                            K_Denominator = [1 -4.48825464011055 7.97295633318904 -6.98919195571153 3.01253347229320 -0.508043209635674 -9.46054293631447e-11 1.81155907741586e-10 -2.09552668075542e-11 -1.38541351445484e-09 6.51896732171504e-09 -2.01631268379705e-08 4.97872366565406e-08 -1.04361724335144e-07 1.90019251788845e-07 -3.02789987135975e-07 4.23547080558947e-07 -5.23107063286633e-07 5.78796240993122e-07 -5.90933953278909e-07 5.83687876272115e-07 -5.87449563212635e-07 6.16624545509239e-07 -6.59780603519169e-07 6.87550063596296e-07 -6.70878296581305e-07 5.97626609977451e-07 -4.78604935047835e-07 3.40971825699377e-07 -2.14212172613813e-07 1.17632031467494e-07 -5.58847655591578e-08 2.26468269180091e-08 -7.63855929526216e-09 2.02186698003906e-09 -3.34459375561136e-10 -2.97501320056686e-11 5.26953894789693e-11 -2.48500884470557e-11 -0.999999999992951 4.48825464010958 -7.97295633318929 6.98919195571278 -3.01253347230011 0.508043209666017];
                        otherwise
                        error('No appropriate filter time constant Tf ={1,3,8} (parDXCPPhaT_cl.Tf)')
                    end
                otherwise
                error('IMC-Controller are only generated for L_S={7,39}')    
            end
        case 'approx_IMC'
            switch B_NumFr
                case 7
                    switch Tf
                        case '1'
                            %Tf=1
                            K_Nominator = [0 6.77529775282132 -24.0822466052188 31.9601217685260 -18.7682622459842 4.11512008912029];
                            K_Denominator = [1 -3.53691622448512 4.66641444028348 -2.72156821256229 0.592069996763923 0];
                        case '3'
                            %Tf=3
                            K_Nominator = [0 0.818877011910587 -2.91063195443529 3.86276883598376 -2.26837241217125 0.497362236342870];
                            K_Denominator = [1 -3.69367096036970 5.10175517665316 -3.12231947609823 0.714235259814764 0];
                        case '8'
                            %Tf=8
                            K_Nominator = [0 0.118251261630870 -0.420314523119963 0.557809392130502 -0.327567993345012 0.0718223995540825];
                            K_Denominator = [1 -3.74546410630640 5.24559559183264 -3.25473122467799 0.754599739151752 0];
                        otherwise
                        error('No appropriate filter time constant Tf ={1,3,8} (parDXCPPhaT_cl.Tf)')
                    end
                case 39
                    switch Tf                 
                        case '1'
                            %Tf=1
                            K_Nominator = [0 32.2831057110573 -125.683935384568 183.474176119244 -119.027424725078 28.9540847526691];
                            K_Denominator = [1 -3.70629666369483 5.13217186751739 -3.14534614264005 0.719470938817486 0];
                        case '3'
                            %Tf=3
                            K_Nominator = [0 3.90180536772069 -15.1904298833245 22.1750822745347 -14.3859097341806 3.49945275762951];
                            K_Denominator = [1 -3.86305139957941 5.59406378989940 -3.59893597321415 0.867923582894160 0];
                        case '8'
                            %Tf=8
                            K_Nominator = [0 0.563446525742075 -2.19359863849069 3.20222868341845 -2.07742060288831 0.505344145199305];
                            K_Denominator = [1 -3.91484454551611 5.74667695088570 -3.74880604991556 0.916973644545979 0];            
                        otherwise
                        error('No appropriate filter time constant Tf ={1,3,8} (parDXCPPhaT_cl.Tf)')
                    end
                otherwise
                error('IMC-Controller are only generated for L_S={7,39}')    
            end

        case 'PIT1'
            switch Tf
                case '1'
                    %Tf=1
                    K_Nominator = [0 1.44352105029607 -1.42908583979311];
                    K_Denominator = [1 -1.75970675828940 0.759706758289376];
                case '3'
                    %Tf=3
                    K_Nominator = [0 0.174467048891578 -0.172722378402663];
                    K_Denominator = [1 -1.91646149417388 0.916461494173877];
                case '8'
                    %Tf=8
                    K_Nominator = [0 0.0251941968627353 -0.0249422548941180];
                    K_Denominator = [1 -1.96825464010938 0.968254640109407];
                otherwise
                error('No appropriate filter time constant Tf ={1,3,8} (parDXCPPhaT_cl.Tf)')
            end
        case 'PI'
            switch Tf
                case '1'
                    %Tf= 1.5835
                    K_Nominator = [7.76509562030564 -7.68744466410254];
                    K_Denominator = [1 -1];
                case '3'
                    %Tf= 4.7506
                    K_Nominator = [2.65841042059604 -2.63182628851887];
                    K_Denominator = [1 -1];
                case '8'
                    %Tf= 12.6683
                    K_Nominator = [1.00530709664882 -0.995254025682334];
                    K_Denominator = [1 -1];
                otherwise
                error('No appropriate filter time constant Tf ={1,3,8} (parDXCPPhaT_cl.Tf)')
            end
        otherwise
            error('No appropriate controller structure (parDXCPPhaT_cl.ControllerStructure)')
    end
    %End controller choice
end

