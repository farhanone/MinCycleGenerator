from minimal_cycles import *
import os
import matplotlib.pyplot as plt
# from mplcursors import cursor, HoverMode

def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def plot_persistence_curves(out_folder_path):
    # l = [6, 11,22, 27,32, 37]
    # l = [6]
    for i in range(1,40):
        # i=11
        f_id = 30000+i
        plot_data = np.loadtxt(os.path.join(out_folder_path, "persistence_plot", pre_keyword + str(f_id) + "_curve.txt"))
        plt.xlim([0.1, 1])
        # plt.legend(str(i))
        plt.plot(plot_data[:,1], plot_data[:,0], label=i)#,color=cmap(i*5a))
        # cursor()
        # plt.show()
        # break

def find_closest(arr, val):
    idx = np.abs(arr - val).argmin()
    # print(arr[idx])
    return idx
def plot_at_thre(val):
    d = []
    t = []
    for i in range(1,42):
        if i==12:
            continue
        f_id = 30000+i
        plot_data = np.loadtxt(os.path.join(out_folder_path, "persistence_plot", pre_keyword + str(f_id) + "_curve.txt"))
        idx = find_closest(plot_data[:, 1], val)
        d.append(plot_data[idx])
        t.append(i)
        # plt.xlim([0, 60])
    #plt.plot(d[], plot_data[:,0], color=cmap(i*6))
        # print(i)
    return np.array(d), t

if __name__ == "__main__":
    
    #### Code to load cuctom data #####
    datatype = "Zero_Link" # "One_Link_61_Degree_ClockWise_Horizontal" #
    # datatype = "One_Link_61_Degree_ClockWise_Horizontal"
    data_folder_path = os.path.join("D:\Projects\GranularMaterials\PhotoElasticityData", datatype ,"FrameWiseResults")
    # out_folder_path = os.path.join("D:\Projects\GranularMaterials\PhotoElasticityData",datatype,"Results")
    out_folder_path = "results"
    pre_keyword = "orientation"
    time_points = 40
    if datatype == "One_Link_61_Degree_ClockWise_Horizontal":
        pre_keyword = "1chain_orientation"
        time_points = 42
     #### Code to load cuctom data #####

    for i in range(1, time_points):
        # i=26
        f_id = 30000+i
        adj_file = os.path.join(data_folder_path, pre_keyword + str(f_id) + "_joAdjacencyAbs.dlm")
        solved_file = os.path.join(data_folder_path, pre_keyword + str(f_id) + "_solved.mat")
        
        
        adj_path = os.path.join(out_folder_path, "adj_matrix")
        create_dir(adj_path)

        adj_loc_radi_to_csv(adj_file, solved_file, os.path.join(adj_path, pre_keyword + str(f_id) + "_adj.csv"))
        df = pd.read_csv(os.path.join(adj_path, pre_keyword + str(f_id) + "_adj.csv"))
        
        VR = os.path.join(out_folder_path, "vietoris_rips")
        create_dir(VR)
        simplex_tree, lookup_id = get_VR(df, os.path.join(VR, pre_keyword + str(f_id) + "_VR.vtp"))

        cycles_dictionary = extract_cycle(simplex_tree, len(df))

        cycles_path = os.path.join(out_folder_path, "cycles")
        create_dir(cycles_path)
        write_cycles_vtp(cycles_dictionary, lookup_id, os.path.join(cycles_path, pre_keyword + str(f_id) + "_cycles.vtp"))

        cycles =[]
        filteration = []
        for key, value in cycles_dictionary.items():
            cycles.append(1)
            filteration.append(value['persistance_value'][0])
        cycles_cum = np.cumsum(cycles)
        
        create_dir(os.path.join(out_folder_path, "persistence_plot"))
        np.savetxt(os.path.join(out_folder_path, "persistence_plot", pre_keyword + str(f_id) + "_curve.txt"), np.stack((cycles_cum, filteration), axis=1))

        # plt.plot(filteration, cycles_cum)
        # plt.savefig(os.path.join(out_folder_path, "persistence_plot", pre_keyword + str(f_id) + "_curve.png"))
        # # plt.clf()


        ########Uncomment########
        ###### Plot Number of Loops vs Filtraion for all time points########
        # plt.figure(figsize=(8, 8))
        # plt.xlabel("Filteration")
        # plt.ylabel("# Loops")
        # plot_persistence_curves(out_folder_path)


        ########Uncomment########
        ###How the number of loops changes w.r.t time while keeping threshold constant.
        #### The threshold is not same where the loops apears for all time points, so we input range of threshold

        # cycle_vs_time, t=plot_at_thre(17)
        # df_plot = pd.DataFrame([cycle_vs_time[:,0].tolist(), cycle_vs_time[:,1].tolist(), t]).transpose()#, columns = [" # Cycles", "Filteration", "Time"])
        # df_plot.columns = ["Cycles", "Filteration", "Time"]

        # plt.figure(figsize=(15, 7))
        # plt.subplot(1, 2, 1)
        # plt.plot(df_plot.Time, df_plot.Cycles)
        # plt.title("Filteration Range {min_: .2f} - {max_: .2f} ".format(min_=min(cycle_vs_time[:,1]), max_=max(cycle_vs_time[:,1])))
        # plt.ylabel("Number of Loops")
        # plt.xlabel("Time")
        # plt.subplot(1, 2, 2)
        # plt.plot(df_plot.Filteration)
        # plt.xlabel("Time")
        # plt.title("Distribution of filteration value")

        # # cursor(hover=True)
        # #cursor()
        # plt.show()