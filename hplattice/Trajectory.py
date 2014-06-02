
class Trajectory(object):
    """docstring for Trajectory"""
    def __init__(self, save_trajectory, trajectory_filename):
        self.save_trajectory = save_trajectory
        self.trajectory_filename = trajectory_filename
        if self.save_trajectory:
            self.output_stream = open(self.trajectory_filename, 'w')
        else:
            self.output_stream = None
        self.frame_num = 0

    def snapshot(self, chain):
        if self.save_trajectory:
            # save chain coords in xyz format
            hp_string = chain.get_hp_string()
            coord_array = chain.get_coord_array()
            coord_strings = \
                ["%s\t%.1f\t%.1f\t%.1f\n" % (hp, row[0], row[1], 0.0) for hp, row in zip(hp_string, coord_array)]
            coord_str = "".join(coord_strings)
            num_atoms = coord_array.shape[0]
            self.output_stream.write("%d\nFrame %d\n" % (num_atoms, self.frame_num))
            self.output_stream.write(coord_str)
            self.frame_num += 1
        else:
            return

    def finalize(self):
        self.output_stream.close()
